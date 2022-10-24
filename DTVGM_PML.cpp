#include <Rcpp.h>
using namespace Rcpp;

/*
 * DTVGM_PML.cpp, C++ function exporting to R
 * Hydrological simulation based on the DTVGM-PML
 * the Distributed Time Variant Gain Model coupled with Penman–Monteith–Leuning (DTVGM-PML)
 *
 * Created by Songzh, 2022-10-24
 *
 * @references
 * [1] Song, Z., Xia, J., Wang, G., She, D., Hu, C., Hong, S., 2022. Regionalization
 *     of hydrological model parameters using gradient boosting machine. Hydrology
 *     and Earth System Sciences 26, 505–524. https://doi.org/10.5194/hess-26-505-2022
 * [2] Zhang, Y., Kong, D., Gan, R., Chiew, F.H.S., McVicar, T.R., Zhang, Q., Yang, Y., 2019.
 *     Coupled estimation of 500 m and 8-day resolution global evapotranspiration and gross
 *     primary production in 2002–2017. Remote Sensing of Environment 222, 165–182.
 *     https://doi.org/10.1016/j.rse.2018.12.031
 * [3] https://github.com/gee-hydro/gee_PML
 */

/* DTVGM_PML
 * data: forcing data including prec, temp, pres, wind, shum, srad, lrad, lai, alb
 * para: model parameters fer0	sls	hc	gsx	cmelt	g1	g2	ks	kr	kg	wm
 * W_init: initial W
 * G_init: initial G
 */
// [[Rcpp::export]]
List DTVGM_PML(NumericMatrix data, NumericVector para, double W_init = 0.0, double G_init = 0.0)
{
  // parameter
  double fer0 = para[0];
  double sls = para[1];
  double hc = para[2];
  double gsx = para[3];
  double cmelt = para[4];
  double g1 = para[5];
  double g2 = para[6];
  double ks = para[7];
  double kr = para[8];
  double kg = para[9];
  double wm = para[10];

  // forcing data length, usually days
  int N = data.nrow();

  // forcing variables
  double prec, temp, pres, wind, shum, srad, lrad, lai, alb;

  // snowmelt variables
  double spk_init = 0.0, rsm_init = 0.0, tt = 0.0, cwh = 0.1, cfr = 0.05;
  double spk0 = spk_init, spk1, rsm0 = rsm_init, rsm1;
  // double pr, ps;
  NumericVector pr(N), smtd(N);

  // interception variables
  // double fveg, sveg, fER, prec_wet, lairef = 5.0;
  NumericVector Ei(N);

  // evapotranspiration variables
  double sigma = 4.9e-9; // Stefan-Boltzmann constant [MJ K-4 m-2 day-1].
  double kmar = 0.4;     // von Karman's constant
  double Zob = 10;       // m, making sure higher than hc, reference height where μ is measured (m).
  double Cp = 0.001013;  // specific heat at constant pressure  [MJ kg-1 K-1] equal to [MJ kg-1 ℃-1]
  double kQ = 0.4488;    // extinction coefficient
  double kA = 0.7;       // the attenuation of net all-wave irradicance, typically about 0.6-0.8 (Denmend, 1976, Kelliher FM et al., (1995))
  double Q50 = 2.6;      // the value of absorbed PAR when gs=gsx/2, [MJ m-2 day-1]
  double D50 = 0.8;      // the humidity deficit at which stomatal conductance is half its maximum value [kpa]
  double lai_ref = 5;    // reference lai
  NumericVector pet(N), Es(N), Ec(N);

  // runoff variables
  NumericVector Rs(N), Rss(N), Rg(N), W(N + 1), G(N + 1);
  W[0] = W_init;
  G[0] = G_init;

  for (int i = 0; i < N; i++)
  {
    // forcing variables
    prec = data(i, 0);
    temp = data(i, 1);
    pres = data(i, 2);
    wind = data(i, 3);
    shum = data(i, 4);
    srad = data(i, 5);
    lrad = data(i, 6);
    lai = data(i, 7);
    alb = data(i, 8);

    double fval_soil = W[i] / wm;

    // snowmelt routinue -------------------------------
    // seperating precipitation
    pr[i] = temp < tt ? 0 : prec;     // rainfall
    double ps = temp < tt ? prec : 0; // snowfall
    if (temp > tt)
    {
      double smt = fmin(spk0 + ps, cmelt * (temp - tt));
      spk1 = spk0 + ps - smt;
      rsm1 = rsm0 + smt;
      double tmp = rsm1 - cwh * spk1;
      if (tmp > 0)
        smtd[i] = tmp;
      else
        smtd[i] = 0;
      rsm1 = rsm1 - smtd[i];
    }
    else
    {
      double rmw = fmin(rsm0, cfr * cmelt * (tt - temp));
      rsm1 = rsm0 - rmw;
      spk1 = spk0 + ps + rmw;
      smtd[i] = 0;
    }
    spk0 = spk1;
    rsm0 = rsm1;
    // snowmelt routinue end ---------------------------

    // interception routine ----------------------------
    double fveg = 1 - exp(-lai / lai_ref);
    double sveg = sls * lai;
    double fER = fveg * fer0;
    double prec_wet = -log(1 - fer0) / fER * sveg;
    if (pr[i] < prec_wet)
      Ei[i] = fveg * pr[i];
    else
      Ei[i] = fveg * prec_wet + fER * (pr[i] - prec_wet);
    Ei[i] = Ei[i] < 0 ? 0 : Ei[i];
    // interception routine end ------------------------

    // actual evapotranspiration routine ---------------
    // Atmosphere parameters
    double lambda = 2.501 - 2.361e-3 * temp;                 // Latent heat of vaporization
    double va = shum * pres / (0.622 + 0.378 * shum);        // actual vapor pressure
    double vs = 0.6108 * exp(17.27 * temp / (temp + 237.3)); // saturated vapour pressure
    double vpd = vs - va;
    double rou_a = 1.293 * 273 / (273 + temp) * pres / 101.3; // density of air [kg/m3]
    double gamma = 0.00163 * pres / lambda;                   // psychrometric constant (S2.9)
    double delta = 4098 * vs / pow(temp + 237.3, 2);          // slope of vapour pressure curve (S2.4)

    // Radiation calculation
    double fv = 1 - exp(-lai / lai_ref);     // vegetation fractional cover
    double se = 0.97 * fv + 0.92 * (1 - fv); // surface emissivity
    double Rns = srad * (1 - alb);
    double Rnl = lrad - se * sigma * pow(temp + 273.2, 4);
    double Rn = Rns + Rnl;
    Rn = Rn < 0 ? 0 : Rn;

    // estimate daily FAO-56 reference crop evapotranspiration
    double u2 = wind * 4.87 / (log(67.8 * Zob - 5.42));
    pet[i] = (0.408 * delta * Rn + gamma * 900 / (temp + 273) * u2 * vpd) / (delta + gamma * (1 + 0.34 * u2));

    // flux density of visible radiation at the top of the canopy (approximately half of incoming solar radiation)
    double Qh = fmax(srad * 0.45, 0);
    double Gc = gsx / kQ * log((Qh + Q50) / (Qh * exp(-kQ * lai) + Q50)) / (1 + delta / D50); // canopy conductance (m/s)
    Gc = fmax(Gc, 1e-6);

    double d = hc * 2 / 3;   // zero plane displacement height (m)
    double zom = 0.123 * hc; // roughness lengths governing transfer
    double zov = 0.1 * zom;  // of momentum and water vapor (m)
    double uz = wind;
    double Ga = kmar * kmar * uz / (log((Zob - d) / zom) * log((Zob - d) / zov)); // aerodynamic conductance (m/s)
    Ga = fmin(Ga, 0.033);

    // Equilibrium evaporation
    // double Eeq = delta / (delta + gamma) * Rn;
    // the fraction of the total available energy (Rn - G) that is absorbed by the soil surface
    double tou = exp(-kA * lai);
    // Transpiration from plant cause by radiation water transfer
    double LEcr = delta / gamma * Rn * (1 - tou) / (delta / gamma + 1 + Ga / Gc);
    // Transpiration from plant cause by aerodynamic water transfer
    double LEca = rou_a * Cp * vpd * Ga / gamma / (delta / gamma + 1 + Ga / Gc);
    LEcr = lai <= 0 ? 0 : LEcr;
    LEca = lai <= 0 ? 0 : LEca;
    double LEc = LEcr + LEca;
    // soil evaporation at equilibrium
    double LEs_eq = delta * Rn * tou / (delta + gamma);
    // evaporation components [mm] [MJ m-2 day-1] / [MJ kg-1] = [kg m-2 day-1]
    double Es_eq = LEs_eq / lambda;
    // double Ecr = LEcr / lambda;
    // double Eca = LEca / lambda;
    Ec[i] = fmax(0, LEc / lambda);
    Es[i] = fmax(0, Es_eq * fval_soil);
    // E[i] = Ei[i] + Ec[i] + Es[i];
    // actual evapotranspiration routine end ------------

    // runoff routine -----------------------------------
    double pe = pr[i] - Ei[i] + smtd[i];
    double et = Es[i] + Ec[i];
    // surface runoff
    double g = g1 * pow(fval_soil, g2);
    g = fmin(g, 1);
    if (pe < et)
      Rs[i] = 0;
    else
      Rs[i] = pe * g;

    double inf = pe - Rs[i];                     // infiltration
    Rss[i] = ks * inf * fval_soil;               // interflow
    Rg[i] = kg * G[i];                           // base flow
    double Gr = kr * fval_soil * (inf - Rss[i]); // groundwater recharge
    double Sr = inf - Rss[i] - Gr;               // soil moisture storage recharge
    W[i + 1] = W[i] + Sr - et;                   // soil moisture storage change
    G[i + 1] = G[i] + Gr - Rg[i];                // groundwater storage change

    if (W[i + 1] < 0)
    {
      G[i + 1] = G[i + 1] + W[i + 1];
      W[i + 1] = 0;
      if (G[i + 1] < 0)
      {
        Es[i] = Es[i] * (1 + G[i + 1] / et);
        Ec[i] = Ec[i] * (1 + G[i + 1] / et);
        G[i + 1] = 0;
      }
    }
    else if (W[i + 1] > wm)
    {
      Rss[i] = Rss[i] + (W[i + 1] - wm);
      W[i + 1] = wm;
    }
    // R[i] = Rs[i] + Rss[i] + Rg[i];
  }

  NumericVector P = pr + smtd;
  NumericVector R = Rs + Rss + Rg;
  NumericVector E = Ei + Es + Ec;

  List out = List::create(Named("P") = P,
                          Named("E") = E,
                          Named("Ei") = Ei,
                          Named("Es") = Es,
                          Named("Ec") = Ec,
                          Named("PET") = pet,
                          Named("R") = R,
                          Named("Rs") = Rs,
                          Named("Rss") = Rss,
                          Named("Rg") = Rg,
                          Named("W") = W,
                          Named("G") = G);
  return (out);
}
