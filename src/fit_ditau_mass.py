import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

ROOT.gROOT.SetBatch(True)


def run_roofit_mass_fit(root_path="plots/ditau_histograms.root"):
    infile = ROOT.TFile.Open(root_path, "READ")
    h_os = infile.Get("h_os")

    m_vis = ROOT.RooRealVar("m_vis", "Visible Di-Tau Mass [GeV]", 30.0, 150.0)
    data_hist = ROOT.RooDataHist("data_hist", "OS Di-Tau Data", ROOT.RooArgList(m_vis), h_os)

    mean = ROOT.RooRealVar("mean", "Peak Mean", 88.0, 70.0, 100.0)
    sigma = ROOT.RooRealVar("sigma", "Peak Width", 6.0, 2.0, 20.0)
    gauss_sig = ROOT.RooGaussian("gauss_sig", "Resonance Peak", m_vis, mean, sigma)

    tau_slope = ROOT.RooRealVar("tau_slope", "Background Slope", -0.02, -0.2, 0.0)
    exp_bkg = ROOT.RooExponential("exp_bkg", "Continuum Background", m_vis, tau_slope)

    n_sig = ROOT.RooRealVar("n_sig", "Signal Yield", 25.0, 0.0, 200.0)
    n_bkg = ROOT.RooRealVar("n_bkg", "Background Yield", 10.0, 0.0, 200.0)
    model = ROOT.RooAddPdf(
        "model", "Signal + Background",
        ROOT.RooArgList(gauss_sig, exp_bkg),
        ROOT.RooArgList(n_sig, n_bkg)
    )

    model.fitTo(data_hist, ROOT.RooFit.PrintLevel(-1))

    x_vals = np.linspace(30.0, 150.0, 300)
    tot_curve = []
    bkg_curve = []
    bin_width = 5.0
    total_expected = n_sig.getVal() + n_bkg.getVal()

    arg_set = ROOT.RooArgSet(m_vis)
    for x in x_vals:
        m_vis.setVal(x)
        tot_pdf = model.getVal(arg_set) * total_expected * bin_width
        bkg_pdf = exp_bkg.getVal(arg_set) * n_bkg.getVal() * bin_width
        tot_curve.append(tot_pdf)
        bkg_curve.append(bkg_pdf)

    bins = np.linspace(30.0, 150.0, 25)
    counts = np.array([h_os.GetBinContent(i) for i in range(1, 25)])

    hep.style.use("CMS")
    fig, ax = plt.subplots(figsize=(8, 6))

    hep.histplot(
        counts, bins, yerr=np.sqrt(np.maximum(counts, 1.0)),
        histtype="errorbar", color="black", label=r"Opposite-Sign Data ($\mu^\pm\tau_h^\mp$)", ax=ax
    )
    ax.plot(x_vals, tot_curve, color="#e42536", lw=2.5, label="RooFit Total Model (S+B)")
    ax.plot(x_vals, bkg_curve, color="#5790fc", linestyle="--", lw=2, label="Exponential Background")

    fit_info = (
        f"Fitted Peak: {mean.getVal():.1f} +/- {mean.getError():.1f} GeV\n"
        f"Resolution: {sigma.getVal():.1f} +/- {sigma.getError():.1f} GeV\n"
        f"Signal Yield: {n_sig.getVal():.1f} +/- {n_sig.getError():.1f} events"
    )
    ax.text(
        0.05, 0.78, fit_info, transform=ax.transAxes, fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc")
    )

    ax.set_xlabel(r"Visible Di-Tau Mass $m_{\mathrm{vis}}(\mu, \tau_h)$ [GeV]")
    ax.set_ylabel("Events / 5 GeV")
    ax.legend(loc="upper right", fontsize=11)
    hep.cms.label("Open Data", data=True, year=2016, ax=ax)

    plt.tight_layout()
    plt.savefig("plots/roofit_mvis_fit.png", dpi=150)
    infile.Close()

    print(f"RooFit Summary -> Mean: {mean.getVal():.2f} GeV | Sigma: {sigma.getVal():.2f} GeV | Signal Yield: {n_sig.getVal():.1f}")


if __name__ == "__main__":
    run_roofit_mass_fit()
