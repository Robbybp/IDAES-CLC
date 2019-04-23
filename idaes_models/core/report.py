"""
Reporting and pretty graphical displays.
"""
__author__ = 'Dan Gunter <dkgunter@lbl.gov>'
__date__ = '5/17/16'

# stdlib
import abc
import math
import pandas as pd
from six import StringIO
# third-party
from IPython.display import display, HTML
from jinja2 import Environment, PackageLoader

class FlowsheetOutput(object):
    """Represents a Flowsheet model's output

    Args:
        model: Completed model
        properties: Model properties
    """
    __metaclass__ = abc.ABCMeta

    SECTIONS = {} # override in derived classes
    PLOT = 'p'
    TABLE = 't'

    def __init__(self, model=None, properties=None):
        """Create with a flowsheet model.

        Args:
            model: pymomo Flowsheet model
            properties: Properties base class
        """
        self._model = model  # TODO [From Qi]: I think this is deprecated code
        assert self._model
        self._prop = properties
        assert self._prop
        self._tables, self._plots = {}, {}

    def report(self, **kw):
        """Show the standard report.

        Args:
            **kw: Items to show. Keywords from `self.SHOW_VALUES`
                  should have string values containing codes from
                  `self.SHOW_*_CODE`

        Returns:
            None
        """
        self._tables, self._plots = self._select_report_sections(kw)
        header = self.header()
        if header:
            display(header)
        for value in self.SECTIONS:
            if value not in self._tables and value not in self._plots:
                continue
            # display title
            title = HTML('<h2>{}</h2>'.format(self.SECTIONS[value]))
            display(title)
            # if there is a table, show it
            if value in self._tables:
                meth = getattr(self, '{}_table'.format(value))
                display(meth(*self._tables[value]))
            # if there is a graph, show it
            if value in self._plots:
                meth = getattr(self, '{}_plot'.format(value))
                plot = meth(*self._plots[value])
                display()
        footer = self.footer()
        if footer:
            display(footer)

    @abc.abstractmethod
    def header(self):
        pass

    @abc.abstractmethod
    def footer(self):
        pass

    def _select_report_sections(self, kw):
        """Get tables and graphs to generate for the
           standard report.
        """
        tables, plots = {}, {}
        for key, value in kw.items():
            key = key.lower()
            if key in self.SECTIONS:
                for code in value:
                    if code == self.PLOT:
                        plots[key] = ()
                    elif code == self.TABLE:
                        tables[key] = ()
                    else:
                        raise ValueError('Unrecognized value code "{}" for {}'
                                         .format(code, key))
            else:
                raise ValueError('Unknown output value "{}"'.format(key))
        return tables, plots

    def style_table(self, df, fmt):
        return df.style.format(fmt)\
                  .set_table_styles([
            dict(selector='td', props=[('border-color', '#CCC'), ('padding', '0 5px 0 5px')]),
            dict(selector='th.row_heading', props=[('color', '#666'), ('font-weight', '300'), ('border', '0')]),
            dict(selector='tr', props=[('border', '0')]),
            dict(selector='', props=[('border', '0')]),
            dict(selector='th.blank', props=[('border', '0')]),
            dict(selector='th.col_heading', props=[('background-color', '#EEE'), ('border-color', '#FFF')])
        ])

class _ignore(object):
    def _repr_html_(self):
        """HTML representation
        """
        env = Environment(loader=PackageLoader('idaes_models.core', 'templates'))
        template = env.get_template('flowsheet_output.html')
        model = self._model
        calculated_values = {
            'co2_capture': 1 - model.absorb.Cu_v[0, 'CO2'].value / model.absorb.Cu_v[1, 'CO2'].value
        }
        params = calculated_values
        params.update(self.show)
        return template.render(model=model, **params)

    def __repr__(self):
        """Default representation, as text.
        """
        R_gas = self._prop.FluidProperties.R_gas
        model = self._model
        ostrm = StringIO()
        ostrm.write("\n\nFLOWS\n")
        ostrm.write("Area: {0} m2".format(model.absorb.A.value))
        ostrm.write("Liquid in: {0} mol/s".format(model.absorb.F_liq1.value))
        ostrm.write("Vapor in: {0} mol/s".format(model.absorb.F_vap1.value))
        ostrm.write("CO2 Capture: {0}".format(
            1 - model.absorb.Cu_v[0, 'CO2'].value
            / model.absorb.Cu_v[1, 'CO2'].value))

        # z-dimension
        ostrm.write("\n\n**z-LOCATION, z is dimnsionless length*\n")
        ostrm.write("   n,        z,      z*L")
        ostrm.write("------------------------")
        for i, z in enumerate(model.absorb.z):
            ostrm.write("{0:4d}, {1:8.4f}, {2:8.4f}"
                  .format(i, z, model.absorb.L.value * z))

        # Temperature profiles.
        ostrm.write("\n\n**TEMPERATURE PROFILES**\n")
        ostrm.write("     T_l,      T_v,   deltaT")
        ostrm.write("-----------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:8.4f}, {1:8.4f}, {2:8.4f}".format(
                model.absorb.T_l[z].value,
                model.absorb.T_v[z].value,
                model.absorb.deltaT[z].value))

        # CO2 profiles.
        ostrm.write("\n\nCO2 PROFILES\n")
        ostrm.write("  x_CO2,   y_CO2,      y*,       y-y*,   alpha,      "
              "C*_LCO2,      He,  He/R/T")
        ostrm.write("-----------------------------------------------------"
              "--------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:7.5f}, {1:7.5f}, {2:10.3e}, {3:7.5f}, {4:7.5f}, "
                  "{5:12.5e}, {6:7.1f}, {7:7.4f}".format(
                model.absorb.prop_l_bulk.y[z, 'CO2'].value,
                model.absorb.prop_v_bulk.y[z, 'CO2'].value,
                model.absorb.y_star[z, 'CO2'].value,
                model.absorb.dy[z, 'CO2'].value,
                model.absorb.prop_l_bulk.alpha[z].value,
                math.exp(model.absorb.logC[z, 'CO2'].value),
                model.absorb.prop_l_bulk.henry[z, 'CO2'].value,
                model.absorb.prop_l_bulk.henry[z, 'CO2'].value
                / R_gas \
                / model.absorb.prop_l_bulk.T[z].value))

        # Equilibrium constants
        ostrm.write("\n\nEQUILIBRIUM CONSTANTS logK, K in appropriate SI units\n")
        ostrm.write("  logK1,   logK2,   logK3,   logK8,   logK9")
        ostrm.write("--------------------------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:7.2f}, {1:7.2f}, {2:7.2f}, {3:7.2f}, {4:7.2f}".format(
                model.absorb.logK[z, 1].value,
                model.absorb.logK[z, 2].value,
                model.absorb.logK[z, 3].value,
                model.absorb.logK[z, 8].value,
                model.absorb.logK[z, 9].value))

        # Chemical Equilibrium
        hs = [s for s in model.absorb.ccomp]
        t = ["{0[" + str(j) + "]:8.2e}" for j, i in enumerate(model.absorb.ccomp)]
        h = ["{:>8}".format(i) for i in hs]
        h = ", ".join(h)
        t = ", ".join(t)

        ostrm.write("\n\nCHEMICAL EQUILIBRIUM CONCENTRATIONS\n")
        ostrm.write(h)
        ostrm.write("---------------------------------------------------"
              "----------------------------")
        for z in model.absorb.z:
            C = [math.exp(model.absorb.logC[z, i].value)
                 for i in model.absorb.ccomp]
            ostrm.write(t.format(C))

        # CHEM EQUILILIBRIUM BALANCES CHECK
        hs = ['sum C', 'sum R', "sum chg"]
        t = ["{0:8.2f}", "{1:8.2f}", "{2:10.2e}"]
        h = ["{0[0]:>8}", "{0[1]:>8}", "{0[2]:>10}"]
        h = ", ".join(h)
        t = ", ".join(t)
        ostrm.write("\n\nCHEMICAL EQUILIBRIUM BALANCES\n")
        ostrm.write("  C = Total Carbon Concentration")
        ostrm.write("  R = Total MEA Species Concenration")
        ostrm.write("  chg = Total Charge Concentration\n")
        ostrm.write(h.format(hs))
        ostrm.write("-------------------------------")
        for z in model.absorb.z:
            sumC = sum(math.exp(model.absorb.logC[z, i].value) for i in
                       ['CO3--', 'CO2', 'HCO3-', 'RNHCOO-'])
            sumR = sum(math.exp(model.absorb.logC[z, i].value) for i in
                       ['RNH2', 'RNH3+', 'RNHCOO-'])
            sumChg = sum(math.exp(model.absorb.logC[z, i].value) for i in
                         ['RNH3+', 'H3O+'])
            sumChg -= sum(math.exp(model.absorb.logC[z, i].value) for i in
                          ['RNHCOO-', 'HCO3-', 'OH-', 'CO3--', 'CO3--'])
            ostrm.write(t.format(sumC, sumR, sumChg))

        # APARENT CONC.
        ostrm.write("\n\nAPPARENT CONCENTRATIONS\n")
        ostrm.write("   L CO2,    L MEA,     V CO2,     V MEA")
        ostrm.write("-----------------------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:8.2f}, {1:8.2f}, {2:9.2e}, {3:9.2e}".format(
                model.absorb.prop_l_bulk.conc_liq[z, 'CO2'].value,
                model.absorb.prop_l_bulk.conc_liq[z, 'MEA'].value,
                model.absorb.prop_v_bulk.conc_vap[z, 'CO2'].value,
                model.absorb.prop_v_bulk.conc_vap[z, 'MEA'].value))

        # total concentration/mole density
        ostrm.write("\n\nTOTAL CONCENTRATION/MOLE DENSITY\n")
        ostrm.write("   Liquid,     Vapor")
        ostrm.write("---------------------")
        for z in model.absorb.z:
            ostrm.write("{0:9.3e}, {1:9.3e}".format(
                model.absorb.prop_l_bulk.conc_liq_tot[z].value,
                model.absorb.prop_v_bulk.conc_vap_tot[z].value))
        # eps profiles.
        ostrm.write("\n\nVOLUME FRACTION AND VELOCITY PROFILES\n")
        ostrm.write(" eps_l,  eps_v,    u_l,    u_v")
        ostrm.write("-------- ----------------------")
        for z in model.absorb.z:
            ostrm.write("{0:5.4f}, {1:5.4f}, {2:5.4f}, {3:5.4f}".format(
                model.absorb.trans.eps_l[z].value,
                model.absorb.trans.eps_v[z].value,
                model.absorb.ul[z].value,
                model.absorb.uv[z].value))

        # Mass/heat transfer
        ostrm.write("\n\nMASS/HEAT LOCAL TRANSFER COEFFICIENTS\n")
        ostrm.write("  k_v_CO2,   k_v_H2O,   k_v_MEA,   k_l_co2,       kap,       koh,          E")
        ostrm.write("-----------------------------------------------------------------------------")
        for i in model.absorb.z:
            ostrm.write("{0:9.3e}, {1:9.3e}, {2:9.3e}, {3:9.3e}, {4:9.3e},  " \
                  "{5:9.3e}, {6:9.3e}".format(
                model.absorb.trans.k_v[i, 'CO2'].value,
                model.absorb.trans.k_v[i, 'H2O'].value,
                model.absorb.trans.k_v[i, 'MEA'].value,
                model.absorb.trans.k_l_CO2[i].value,
                model.absorb.trans.kap[i].value,
                model.absorb.trans.koh[i].value,
                model.absorb.trans.E[i].value))

        # Mass/heat transfer
        ostrm.write("\n\nMASS/HEAT OVERALL TRANSFER COEFFICIENTS\n")
        ostrm.write("  K_v_CO2,   K_v_H2O,   K_v_MEA,         h,       ae/a")
        ostrm.write("-------------------------------------------------------")
        for i in model.absorb.z:
            ostrm.write("{0:9.3e}, {1:9.3e}, {2:9.3e}, {3:9.3e},  " \
                  "{4:9.3e}".format(
                model.absorb.trans.K_v[i, 'CO2'].value,
                model.absorb.trans.K_v[i, 'H2O'].value,
                model.absorb.trans.K_v[i, 'MEA'].value,
                model.absorb.trans.h_vl[i].value,
                model.absorb.trans.ae[i].value / model.absorb.trans.a.value))

        # Mass/heat transfer
        ostrm.write("\n\nMASS AND HEAT TRANSFER\n")
        ostrm.write("      N_CO2,      N_H2O,      N_MEA,        Ntotal,           Q")
        ostrm.write("----------------------------------------------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:11.3e}, {1:11.3e}, {2:11.3e}, {3:11.3e}, {4:11.3e}".format(
                model.absorb.N_v[z, 'CO2'].value,
                model.absorb.N_v[z, 'H2O'].value,
                model.absorb.N_v[z, 'MEA'].value,
                model.absorb.N_v[z, 'CO2'].value \
                + model.absorb.N_v[z, 'H2O'].value \
                + model.absorb.N_v[z, 'MEA'].value,
                model.absorb.Q_v[z].value))

        # Viscosity and Diffusivity
        ostrm.write("\n\nDIFFUSIVITY, VISCOSITY, DENSITY, AND SURFACE TENSION\n")
        ostrm.write("    D_l_CO2,    D_v_CO2,         mu_l,        mu_v,       rho_liq,       rho_vap,       sigma")
        ostrm.write("------------------------------------------------------------------------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:11.3e}, {1:11.3e}, {2:11.3e}, {3:11.3e}" \
                  ", {4:11.3e}, {5:11.3e}, {6:11.3e}".format(
                model.absorb.prop_l_bulk.D_liq[z, 'CO2'].value,
                model.absorb.prop_v_bulk.D_vap[z, 'CO2'].value,
                model.absorb.prop_l_bulk.mu_liq[z].value,
                model.absorb.prop_v_bulk.mu_vap[z].value,
                model.absorb.prop_l_bulk.rho_liq[z].value,
                model.absorb.prop_v_bulk.rho_vap[z].value,
                model.absorb.prop_l_bulk.sigma[z].value))

        #
        ostrm.write("\n\nLIQUID FLOWS (mol/m^2/s)\n")
        ostrm.write("  Cu_l_CO2,   Cu_l_H2O,   Cu_l_MEA")
        ostrm.write("-----------------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:10.3e}, {1:10.3e}, {2:10.3e}".format(
                model.absorb.Cu_l[z, 'CO2'].value,
                model.absorb.Cu_l[z, 'H2O'].value,
                model.absorb.Cu_l[z, 'MEA'].value))

        #
        ostrm.write("\n\nVAPOR FLOWS (mol/m^2/s)\n")
        ostrm.write("  Cu_v_CO2,   Cu_v_H2O,   Cu_v_MEA,    Cu_v_N2")
        ostrm.write("-----------------------------------------------")
        for z in model.absorb.z:
            ostrm.write("{0:10.3e}, {1:10.3e}, {2:10.3e}, {3:10.3e}".format(
                model.absorb.Cu_v[z, 'CO2'].value,
                model.absorb.Cu_v[z, 'H2O'].value,
                model.absorb.Cu_v[z, 'MEA'].value,
                model.absorb.Cu_v[z, 'N2'].value))
        return ostrm.getvalue()

