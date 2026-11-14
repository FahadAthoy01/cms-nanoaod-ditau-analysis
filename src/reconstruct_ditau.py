import ROOT
import yaml
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

ROOT.gInterpreter.Declare("""
#include <cmath>
#include "ROOT/RVec.hxx"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "TVector2.h"

using namespace ROOT::VecOps;

float compute_mt(float mu_pt, float mu_phi, float met_pt, float met_phi) {
    float dphi = TVector2::Phi_mpi_pi(mu_phi - met_phi);
    return std::sqrt(2.0f * mu_pt * met_pt * (1.0f - std::cos(dphi)));
}

float compute_delta_r(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    return std::sqrt(deta * deta + dphi * dphi);
}

float compute_vis_mass(float pt1, float eta1, float phi1, float m1,
                       float pt2, float eta2, float phi2, float m2) {
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, m1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, m2);
    return (p1 + p2).M();
}
""")


def run_ditau_selection(config_path="configs/workflow.yaml"):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    url = cfg["dataset"]["xrootd_url"]
    max_ev = cfg["dataset"]["max_events"]
    sel = cfg["selection"]

    # Range() requires single-threaded RDataFrame execution
    df = ROOT.RDataFrame("Events", url).Range(max_ev)

    df_pairs = (
        df.Filter("nMuon >= 1 && nTau >= 1", "Pre-selection: >=1 Muon and >=1 Tau")
          .Define(
              "mu_mask",
              f"Muon_pt > {sel['muon_pt_min']} && abs(Muon_eta) < {sel['muon_eta_max']} "
              f"&& Muon_tightId == true && Muon_pfRelIso04_all < {sel['muon_iso_max']}"
          )
          .Define(
              "tau_mask",
              f"Tau_pt > {sel['tau_pt_min']} && abs(Tau_eta) < {sel['tau_eta_max']}"
          )
          .Filter("Sum(mu_mask) >= 1 && Sum(tau_mask) >= 1", "Object kinematics and isolation")
          .Define("mu_idx", "ArgMax(Muon_pt * mu_mask)")
          .Define("tau_idx", "ArgMax(Tau_pt * tau_mask)")
          .Define("mu_pt", "Muon_pt[mu_idx]")
          .Define("mu_eta", "Muon_eta[mu_idx]")
          .Define("mu_phi", "Muon_phi[mu_idx]")
          .Define("mu_mass", "Muon_mass[mu_idx]")
          .Define("mu_q", "Muon_charge[mu_idx]")
          .Define("tau_pt", "Tau_pt[tau_idx]")
          .Define("tau_eta", "Tau_eta[tau_idx]")
          .Define("tau_phi", "Tau_phi[tau_idx]")
          .Define("tau_mass", "Tau_mass[tau_idx]")
          .Define("tau_q", "Tau_charge[tau_idx]")
          .Define("dR", "compute_delta_r(mu_eta, mu_phi, tau_eta, tau_phi)")
          .Filter(f"dR > {sel['delta_r_min']}", "DeltaR separation > 0.5")
    )

    df_veto = (
        df_pairs.Define(
            "n_bjet",
            f"Sum(Jet_pt > 20.0 && abs(Jet_eta) < 2.4 && Jet_btagCSVV2 > {sel['btag_csvv2_wp']})"
        )
        .Filter("n_bjet == 0", "b-jet veto for ttbar suppression")
        .Define("mt_mu_met", "compute_mt(mu_pt, mu_phi, MET_pt, MET_phi)")
        .Define("m_vis", "compute_vis_mass(mu_pt, mu_eta, mu_phi, mu_mass, tau_pt, tau_eta, tau_phi, tau_mass)")
    )

    df_sr = df_veto.Filter(f"mt_mu_met < {sel['mt_max']}", "Low mT signal region")
    df_os = df_sr.Filter("mu_q * tau_q < 0", "Opposite-sign signal region")
    df_ss = df_sr.Filter("mu_q * tau_q > 0", "Same-sign QCD control region")

    h_os_ptr = df_os.Histo1D(("h_os", "OS Visible Mass", 24, 30, 150), "m_vis")
    h_ss_ptr = df_ss.Histo1D(("h_ss", "SS Visible Mass", 24, 30, 150), "m_vis")
    h_mt_ptr = df_veto.Histo1D(("h_mt", "Transverse Mass", 30, 0, 150), "mt_mu_met")
    cutflow = df_sr.Report()

    h_os = h_os_ptr.GetValue()
    h_ss = h_ss_ptr.GetValue()
    h_mt = h_mt_ptr.GetValue()

    print("\n--- CMS NanoAOD Cutflow Summary ---")
    cutflow.Print()

    out_root = ROOT.TFile("plots/ditau_histograms.root", "RECREATE")
    h_os.Write()
    h_ss.Write()
    h_mt.Write()
    out_root.Close()

    hep.style.use("CMS")
    bins_m = np.linspace(30, 150, 25)
    os_counts = np.array([h_os.GetBinContent(i) for i in range(1, 25)])
    ss_counts = np.array([h_ss.GetBinContent(i) for i in range(1, 25)])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    hep.histplot(
        ss_counts, bins_m, histtype="fill", color="#f89c20", alpha=0.75,
        label=r"Same-Sign QCD Estimate ($\mu^\pm\tau_h^\pm$)", ax=ax1
    )
    hep.histplot(
        os_counts, bins_m, yerr=np.sqrt(np.maximum(os_counts, 1.0)),
        histtype="errorbar", color="black", label=r"Opposite-Sign Data ($\mu^\pm\tau_h^\mp$)", ax=ax1
    )
    ax1.set_xlabel(r"Visible Di-Tau Mass $m_{\mathrm{vis}}(\mu, \tau_h)$ [GeV]")
    ax1.set_ylabel("Events / 5 GeV")
    ax1.legend(loc="upper right", fontsize=12)
    hep.cms.label("Open Data", data=True, year=2016, ax=ax1)

    bins_mt = np.linspace(0, 150, 31)
    mt_counts = np.array([h_mt.GetBinContent(i) for i in range(1, 31)])
    hep.histplot(
        mt_counts, bins_mt, histtype="step", lw=2, color="#5790fc",
        label=r"Selected $\mu + \tau_h$ Pairs", ax=ax2
    )
    ax2.axvline(sel["mt_max"], color="#e42536", linestyle="--", lw=2, label=r"$m_T < 40$ GeV Cut")
    ax2.set_xlabel(r"Transverse Mass $m_T(\mu, p_T^{\mathrm{miss}})$ [GeV]")
    ax2.set_ylabel("Events / 5 GeV")
    ax2.legend(loc="upper right", fontsize=12)
    hep.cms.label("Open Data", data=True, year=2016, ax=ax2)

    plt.tight_layout()
    plt.savefig("plots/cms_ditau_reconstruction.png", dpi=150)

    os_total = int(np.sum(os_counts))
    ss_total = int(np.sum(ss_counts))
    print(f"\nSignal Region Yield -> OS Events: {os_total} | SS QCD Events: {ss_total} | Net OS-SS Excess: {os_total - ss_total}")


if __name__ == "__main__":
    run_ditau_selection()
