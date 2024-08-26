#include "WireCellSigProc/RefactorSigProc.h"
#include "ROI_formation.h"
#include "ROI_refinement.h"


#include "WireCellAux/DftTools.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IFilterWaveform.h"
#include "WireCellIface/IChannelResponse.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/FFTBestLength.h"
#include "WireCellUtil/Waveform.h"

#include "WireCellUtil/NamedFactory.h"


WIRECELL_FACTORY(RefactorSigProc, WireCell::SigProc::RefactorSigProc,
                 WireCell::INamed,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SigProc;

using WireCell::Aux::DftTools::fwd;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv;
using WireCell::Aux::DftTools::inv_c2r;

RefactorSigProc::~RefactorSigProc() {}


RefactorSigProc::RefactorSigProc(    )
  : Aux::Logger("RefactorSigProc", "sigproc")
{
    // get wires for each plane

    // std::cout << m_anode->channels().size() << " " << nwire_u << " " << nwire_v << " " << nwire_w << std::endl;
}
void RefactorSigProc::load_data(const input_pointer& in, int plane)
{
    m_r_data[plane] = Array::array_xxf::Zero(m_fft_nwires[plane], m_fft_nticks);

    auto traces = in->traces();

    auto& bad = m_wanmm["bad"];
    int nbad = 0;

    for (auto trace : *traces.get()) {
        int wct_channel_ident = trace->channel();
        OspChan och = m_channel_map[wct_channel_ident];
        if (plane != och.plane) {
            continue;  // we'll catch it in another call to load_data
        }

        // fixme: this code uses tbin() but other places in this file will barf if tbin!=0.
        int tbin = trace->tbin();
        auto const& charges = trace->charge();
        const int ntbins = std::min((int) charges.size(), m_nticks);
        for (int qind = 0; qind < ntbins; ++qind) {
            const float q = charges[qind];
            m_r_data[plane](och.wire + m_pad_nwires[plane], tbin + qind) = q;
        }

        // ensure dead channels are indeed dead ...
        auto const& badch = bad.find(och.channel);
        if (badch == bad.end()) {
            continue;
        }

        auto const& binranges = badch->second;
        for (auto const& br : binranges) {
            ++nbad;
            for (int i = br.first; i != br.second; ++i) {
                m_r_data[plane](och.wire + m_pad_nwires[plane], i) = 0;
            }
        }
    }
    //rebase for this plane
    /*if (std::find(m_rebase_planes.begin(), m_rebase_planes.end(), plane) != m_rebase_planes.end()) {
        log->debug("rebase_waveform for plane {} with m_rebase_nbins = {}", plane, m_rebase_nbins);
        rebase_waveform(m_r_data[plane],m_rebase_nbins);
    }*/

    log->debug("call={} load plane index: {}, ntraces={}, input bad regions: {}",
               m_count, plane, traces->size(), nbad);
    //check_data(plane, "load data");
}
void RefactorSigProc::configure(const WireCell::Configuration& config)
{
    m_sparse = get(config, "sparse", false);
    m_load_fr_with_plane_ident = get(config, "load_fr_with_plane_ident", false);

    m_fine_time_offset = get(config, "ftoffset", m_fine_time_offset);
    m_coarse_time_offset = get(config, "ctoffset", m_coarse_time_offset);
    m_anode_tn = get(config, "anode", m_anode_tn);

    std::string dft_tn = get<std::string>(config, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
    m_verbose = get(config, "verbose", 0);

    // m_nticks = get(config,"nticks",m_nticks);
    if (!config["nticks"].isNull()) {
        log->warn("config: no setting \"nticks\", ignoring value {}", config["nticks"].asInt());
    }
    // m_period = get(config,"period",m_period);
    if (!config["period"].isNull()) {
        log->warn("config: no setting \"period\", ignoring value {}", config["period"].asDouble());
    }

    m_fft_flag = get(config, "fft_flag", m_fft_flag);
    if (m_fft_flag) {
      m_fft_flag = 0;
      log->warn("config: fft_flag option is broken, will use native array sizes");
    }
    m_elecresponse_tn = get(config, "elecresponse", m_elecresponse_tn);
    m_gain = get(config, "gain", m_gain);
    m_shaping_time = get(config, "shaping", m_shaping_time);
    m_inter_gain = get(config, "postgain", m_inter_gain);
    m_ADC_mV = get(config, "ADC_mV", m_ADC_mV);

    m_per_chan_resp = get(config, "per_chan_resp", m_per_chan_resp);
    m_field_response = get(config, "field_response", m_field_response);

    m_th_factor_ind = get(config, "troi_ind_th_factor", m_th_factor_ind);
    m_th_factor_col = get(config, "troi_col_th_factor", m_th_factor_col);
    m_pad = get(config, "troi_pad", m_pad);
    m_asy = get(config, "troi_asy", m_asy);
    m_rebin = get(config, "lroi_rebin", m_rebin);
    m_l_factor = get(config, "lroi_th_factor", m_l_factor);
    m_l_max_th = get(config, "lroi_max_th", m_l_max_th);
    m_l_factor1 = get(config, "lroi_th_factor1", m_l_factor1);
    m_l_short_length = get(config, "lroi_short_length", m_l_short_length);
    m_l_jump_one_bin = get(config, "lroi_jump_one_bin", m_l_jump_one_bin);

    m_r_th_factor = get(config, "r_th_factor", m_r_th_factor);
    m_r_fake_signal_low_th = get(config, "r_fake_signal_low_th", m_r_fake_signal_low_th);
    m_r_fake_signal_high_th = get(config, "r_fake_signal_high_th", m_r_fake_signal_high_th);
    m_r_fake_signal_low_th_ind_factor =
        get(config, "r_fake_signal_low_th_ind_factor", m_r_fake_signal_low_th_ind_factor);
    m_r_fake_signal_high_th_ind_factor =
        get(config, "r_fake_signal_high_th_ind_factor", m_r_fake_signal_high_th_ind_factor);
    m_r_pad = get(config, "r_pad", m_r_pad);
    m_r_break_roi_loop = get(config, "r_break_roi_loop", m_r_break_roi_loop);
    m_r_th_peak = get(config, "r_th_peak", m_r_th_peak);
    m_r_sep_peak = get(config, "r_sep_peak", m_r_sep_peak);
    m_r_low_peak_sep_threshold_pre = get(config, "r_low_peak_sep_threshold_pre", m_r_low_peak_sep_threshold_pre);
    m_r_max_npeaks = get(config, "r_max_npeaks", m_r_max_npeaks);
    m_r_sigma = get(config, "r_sigma", m_r_sigma);
    m_r_th_percent = get(config, "r_th_percent", m_r_th_percent);

    if (config.isMember("process_planes")) {
        m_process_planes.clear();
        for (auto jplane : config["process_planes"]) {
            m_process_planes.push_back(jplane.asInt());
        }
    }

    m_charge_ch_offset = get(config, "charge_ch_offset", m_charge_ch_offset);

    m_wiener_tag = get(config, "wiener_tag", m_wiener_tag);
    // m_wiener_threshold_tag = get(config, "wiener_threshold_tag", m_wiener_threshold_tag);
    if (! config["wiener_threshold_tag"].isNull()) {
        log->warn("The 'wiener_threshold_tag' is obsolete, thresholds in summary on 'wiener' tagged traces");
    }
    m_decon_charge_tag = get(config, "decon_charge_tag", m_decon_charge_tag);
    m_gauss_tag = get(config, "gauss_tag", m_gauss_tag);
    m_frame_tag = get(config, "frame_tag", m_frame_tag);

    m_use_roi_debug_mode = get(config, "use_roi_debug_mode", m_use_roi_debug_mode);
    m_use_roi_refinement = get(config, "use_roi_refinement", m_use_roi_refinement);
    m_tight_lf_tag = get(config, "tight_lf_tag", m_tight_lf_tag);
    m_loose_lf_tag = get(config, "loose_lf_tag", m_loose_lf_tag);
    m_cleanup_roi_tag = get(config, "cleanup_roi_tag", m_cleanup_roi_tag);
    m_break_roi_loop1_tag = get(config, "break_roi_loop1_tag", m_break_roi_loop1_tag);
    m_break_roi_loop2_tag = get(config, "break_roi_loop2_tag", m_break_roi_loop2_tag);
    m_shrink_roi_tag = get(config, "shrink_roi_tag", m_shrink_roi_tag);
    m_extend_roi_tag = get(config, "extend_roi_tag", m_extend_roi_tag);

    m_use_multi_plane_protection = get<bool>(config, "use_multi_plane_protection", m_use_multi_plane_protection);
    m_mp3_roi_tag = get(config, "mp3_roi_tag", m_mp3_roi_tag);
    m_mp2_roi_tag = get(config, "mp2_roi_tag", m_mp2_roi_tag);
    m_mp_th1 = get(config, "mp_th1", m_mp_th1);
    m_mp_th2 = get(config, "mp_th2", m_mp_th2);
    m_mp_tick_resolution = get(config, "mp_tick_resolution", m_mp_tick_resolution);
    

    if (config.isMember("rebase_planes")) {
       m_rebase_planes.clear();
       for (auto jplane : config["rebase_planes"]) {
           m_rebase_planes.push_back(jplane.asInt());
       }
    }
    m_rebase_nbins = get(config, "rebase_nbins", m_rebase_nbins);

    m_isWrapped = get<bool>(config, "isWrapped", m_isWrapped);

    // this throws if not found
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

    //
    m_elecresponse = Factory::find_tn<IWaveform>(m_elecresponse_tn);

    // Build up the channel map.  The OSP channel must run contiguously
    // first up the U, then V, then W "wires".  Ie, face-major order,
    // but we have plane-major order so make a temporary collection.
    IChannel::vector plane_channels[3];
    std::stringstream ss;
    ss << "config: internal channel map for tags: gauss:\"" << m_gauss_tag << "\", wiener:\"" << m_wiener_tag
       << "\", frame:\"" << m_frame_tag << "\"\n";

    // fixme: this loop is now available as Aux::plane_channels()
    for (auto face : m_anode->faces()) {
        if (!face) {   // A null face means one sided AnodePlane.
            continue;  // Can be "back" or "front" face.
        }
        for (auto plane : face->planes()) {
            int plane_index = plane->planeid().index();
            auto& pchans = plane_channels[plane_index];
            // These IChannel vectors are ordered in same order as wire-in-plane.
            const auto& ichans = plane->channels();
            // Append
            pchans.reserve(pchans.size() + ichans.size());
            pchans.insert(pchans.end(), ichans.begin(), ichans.end());
            ss << "\tpind" << plane_index << " "
               << "aid" << m_anode->ident() << " "
               << "fid" << face->ident() << " "
               << "pid" << plane->ident() << " "
               << "cid" << ichans.front()->ident() << " -> cid" << ichans.back()->ident() << ", "
               << "cind" << ichans.front()->index() << " -> cind" << ichans.back()->index() << ", "
               << "(n=" << pchans.size() << ")\n";
        }
    }
    log->debug(ss.str());

    int osp_channel_number = 0;
    for (int iplane = 0; iplane < 3; ++iplane) {
        m_nwires[iplane] = plane_channels[iplane].size();
        int osp_wire_number = 0;
        // note the order here is the IChannel::index or Wire Attachment Number
        for (auto ichan : plane_channels[iplane]) {
            const int wct_chan_ident = ichan->ident();
            OspChan och(osp_channel_number, osp_wire_number, iplane, wct_chan_ident);
            // std::cout << "[hyu1]chmap: " << wct_chan_ident << " " << iplane << " " << osp_channel_number << " " <<
            // osp_wire_number << std::endl;
            m_roi_ch_ch_ident[osp_channel_number] = wct_chan_ident;
            m_channel_map[wct_chan_ident] = och;     // we could save some space by storing
            m_channel_range[iplane].push_back(och);  // wct ident here instead of a whole och.
            ++osp_wire_number;
            ++osp_channel_number;
        }
    }
}
WireCell::Configuration RefactorSigProc::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use
    cfg["verbose"] = 0;         // larger is more more logging 
    cfg["ftoffset"] = m_fine_time_offset;
    cfg["ctoffset"] = m_coarse_time_offset;
    // cfg["nticks"] = m_nticks;
    // cfg["period"] = m_period;

    //cfg["fft_flag"] = m_fft_flag;
    cfg["fft_flag"] = 0;

    cfg["elecresponse"] = m_elecresponse_tn;
    cfg["gain"] = m_gain;
    cfg["shaping"] = m_shaping_time;
    cfg["inter_gain"] = m_inter_gain;
    cfg["ADC_mV"] = m_ADC_mV;

    cfg["per_chan_resp"] = m_per_chan_resp;
    cfg["field_response"] = m_field_response;

    cfg["troi_ind_th_factor"] = m_th_factor_ind;
    cfg["troi_col_th_factor"] = m_th_factor_col;
    cfg["troi_pad"] = m_pad;
    cfg["troi_asy"] = m_asy;
    cfg["lroi_rebin"] = m_rebin;
    cfg["lroi_th_factor"] = m_l_factor;
    cfg["lroi_max_th"] = m_l_max_th;
    cfg["lroi_th_factor1"] = m_l_factor1;
    cfg["lroi_short_length"] = m_l_short_length;
    cfg["lroi_jump_one_bin"] = m_l_jump_one_bin;

    cfg["r_th_factor"] = m_r_th_factor;
    cfg["r_fake_signal_low_th"] = m_r_fake_signal_low_th;
    cfg["r_fake_signal_high_th"] = m_r_fake_signal_high_th;
    cfg["r_fake_signal_low_th_ind_factor"] = m_r_fake_signal_low_th_ind_factor;
    cfg["r_fake_signal_high_th_ind_factor"] = m_r_fake_signal_high_th_ind_factor;
    cfg["r_pad"] = m_r_pad;
    cfg["r_break_roi_loop"] = m_r_break_roi_loop;
    cfg["r_th_peak"] = m_r_th_peak;
    cfg["r_sep_peak"] = m_r_sep_peak;
    cfg["r_low_peak_sep_threshold_pre"] = m_r_low_peak_sep_threshold_pre;
    cfg["r_max_npeaks"] = m_r_max_npeaks;
    cfg["r_sigma"] = m_r_sigma;
    cfg["r_th_precent"] = m_r_th_percent;

    // cfg["process_planes"] = Json::arrayValue;

    // fixme: unused?
    cfg["charge_ch_offset"] = m_charge_ch_offset;

    cfg["wiener_tag"] = m_wiener_tag;
    // cfg["wiener_threshold_tag"] = m_wiener_threshold_tag;
    cfg["decon_charge_tag"] = m_decon_charge_tag;
    cfg["gauss_tag"] = m_gauss_tag;
    cfg["frame_tag"] = m_frame_tag;

    cfg["use_roi_debug_mode"] = m_use_roi_debug_mode;  // default false
    cfg["use_roi_refinement"] = m_use_roi_refinement;  // default true
    cfg["tight_lf_tag"] = m_tight_lf_tag;
    cfg["loose_lf_tag"] = m_loose_lf_tag;
    cfg["cleanup_roi_tag"] = m_cleanup_roi_tag;
    cfg["break_roi_loop1_tag"] = m_break_roi_loop1_tag;
    cfg["break_roi_loop2_tag"] = m_break_roi_loop2_tag;
    cfg["shrink_roi_tag"] = m_shrink_roi_tag;
    cfg["extend_roi_tag"] = m_extend_roi_tag;

    cfg["use_multi_plane_protection"] = m_use_multi_plane_protection;  // default false
    cfg["mp3_roi_tag"] = m_mp3_roi_tag;
    cfg["mp2_roi_tag"] = m_mp2_roi_tag;
    cfg["mp_th1"] = m_mp_th1;
    cfg["mp_th2"] = m_mp_th2;
    cfg["mp_tick_resolution"] = m_mp_tick_resolution;
    
    cfg["rebase_nbins"] = m_rebase_nbins;    

    cfg["isWarped"] = m_isWrapped;  // default false

    cfg["sparse"] = false;
    cfg["load_fr_with_plane_ident"] = false;

    return cfg;
}


bool RefactorSigProc::operator()(const input_pointer& in, output_pointer& out) {
  ITrace::vector* itraces = new ITrace::vector;  // will become shared_ptr.
  auto sframe = new Aux::SimpleFrame(in->ident(), in->time(), ITrace::shared_vector(itraces), in->tick(), in->masks());
  sframe->tag_frame(m_frame_tag);
  out = IFrame::pointer(sframe);
  return true;
}
