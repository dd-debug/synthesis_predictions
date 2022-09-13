'''
Created on Jan 20, 2021

@author: jiadongc
'''

from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
import plotly.graph_objects as go
import json
from pymatgen.analysis.reaction_calculator import ComputedReaction
with open("C:/Users/jdche/anaconda3/envs/my_pymatgen/lib/site-packages/pymatgen/util/plotly_pd_layouts.json", "r") as f:
    plotly_layouts = json.load(f)


class Inter_PDPlotter(PDPlotter):
    def __init__(
        self,
        phasediagram: PhaseDiagram,
#         new_entries,
        reactions_dict,
        target_entry = None,
        show_unstable: float = 0.2,
        backend: str = "plotly",
        **plotkwargs,
    ):
        super().__init__(phasediagram,
        show_unstable = show_unstable,
        backend = backend,
        **plotkwargs)
        self.reactions_dict = reactions_dict
        self.target_entry = target_entry
        
    def get_x_y_values(self):
        x, y, z, text, textpositions = [], [], [], [], []
        stable_labels_plot = None
        min_energy_x = None
        offset_2d = 0.005  # extra distance to offset label position for clarity
        offset_3d = 0.01

        energy_offset = -0.1 * self._min_energy
#         for a,e in list(self.pd_plot_data[1].items()):
#             print("haha",e.name)
#         for e in self.new_entries:
#             print("gaga",e.name)
        if self._dim == 2:
            min_energy_x = min(list(self.pd_plot_data[1].keys()), key=lambda c: c[1])[0]
        
        for coords, entry in self.pd_plot_data[1].items():
            if entry.composition.is_element:  # taken care of by other function
                continue
            x_coord = coords[0]
            y_coord = coords[1]
            textposition = None

            if self._dim == 2:
                textposition = "bottom left"
                if x_coord >= min_energy_x:
                    textposition = "bottom right"
                    x_coord += offset_2d
                else:
                    x_coord -= offset_2d
                y_coord -= offset_2d
            elif self._dim == 3:
                textposition = "middle center"
                if coords[0] > 0.5:
                    x_coord += offset_3d
                else:
                    x_coord -= offset_3d
                if coords[1] > 0.866 / 2:
                    y_coord -= offset_3d
                else:
                    y_coord += offset_3d

                z.append(self._pd.get_form_energy_per_atom(entry) + energy_offset)

            elif self._dim == 4:
                x_coord = x_coord - offset_3d
                y_coord = y_coord - offset_3d
                textposition = "bottom right"
                z.append(coords[2])

            x.append(x_coord)
            y.append(y_coord)
            textpositions.append(textposition)

#             comp = entry.composition
#             if hasattr(entry, "original_entry"):
#                 comp = entry.original_entry.composition
# 
#             formula = list(comp.reduced_formula)
#             text.append(self._htmlize_formula(formula))
            if self.target_entry != None:
                if entry.name == self.target_entry.name:
                    text.append(self.target_entry.name)
                else:
                    text.append("")
            else:

                if self.reactions_dict[entry.name]._reactant_entries == \
                self.reactions_dict[entry.name]._product_entries:
                    text.append(entry.name)
                else: 
                    ''' If the entry is made manually, do not text'''
                    text.append("")
#                     tx = self.reactions_dict[entry.name].__str__().split("->")[-1]
#                     text.append(tx)
        print(x)
        print(y)  
        return(x,y)

    def _create_plotly_stable_labels(self, label_stable=True):
        """
        Creates a (hidable) scatter trace containing labels of stable phases.
        Contains some functionality for creating sensible label positions.

        :return: go.Scatter (or go.Scatter3d) plot
        """
        x, y, z, text, textpositions = [], [], [], [], []
        stable_labels_plot = None
        min_energy_x = None
        offset_2d = 0.005  # extra distance to offset label position for clarity
        offset_3d = 0.01

        energy_offset = -0.1 * self._min_energy
#         for a,e in list(self.pd_plot_data[1].items()):
#             print("haha",e.name)
#         for e in self.new_entries:
#             print("gaga",e.name)
        if self._dim == 2:
            min_energy_x = min(list(self.pd_plot_data[1].keys()), key=lambda c: c[1])[0]
        
        for coords, entry in self.pd_plot_data[1].items():
            if entry.composition.is_element:  # taken care of by other function
                continue
            x_coord = coords[0]
            y_coord = coords[1]
            textposition = None

            if self._dim == 2:
                textposition = "bottom left"
                if x_coord >= min_energy_x:
                    textposition = "bottom right"
                    x_coord += offset_2d
                else:
                    x_coord -= offset_2d
                y_coord -= offset_2d
            elif self._dim == 3:
                textposition = "middle center"
                if coords[0] > 0.5:
                    x_coord += offset_3d
                else:
                    x_coord -= offset_3d
                if coords[1] > 0.866 / 2:
                    y_coord -= offset_3d
                else:
                    y_coord += offset_3d

                z.append(self._pd.get_form_energy_per_atom(entry) + energy_offset)

            elif self._dim == 4:
                x_coord = x_coord - offset_3d
                y_coord = y_coord - offset_3d
                textposition = "bottom right"
                z.append(coords[2])

            x.append(x_coord)
            y.append(y_coord)
            textpositions.append(textposition)

#             comp = entry.composition
#             if hasattr(entry, "original_entry"):
#                 comp = entry.original_entry.composition
# 
#             formula = list(comp.reduced_formula)
#             text.append(self._htmlize_formula(formula))
            if self.target_entry is not None:
                if entry.name == self.target_entry.name:
                    text.append(self.target_entry.name)
                else:
                    text.append("")
            else:

                if self.reactions_dict[entry.name]._reactant_entries == \
                self.reactions_dict[entry.name]._product_entries:
                    text.append(entry.name)
                else: 
                    ''' If the entry is made manually, do not text'''
                    text.append("")
#                     tx = self.reactions_dict[entry.name].__str__().split("->")[-1]
#                     text.append(tx)
        print(x)
        print(y)             
        visible = True
        if not label_stable or self._dim == 4:
            visible = "legendonly"

        plot_args = dict(
            text=text,
            textposition=textpositions,
            mode="text",
            name="Labels (stable)",
            hoverinfo="skip",
            opacity=1.0,
            visible=visible,
            showlegend=True,
        )

        if self._dim == 2:
            stable_labels_plot = go.Scatter(x=x, y=y, **plot_args)
        elif self._dim == 3:
            stable_labels_plot = go.Scatter3d(x=y, y=x, z=z, **plot_args)
        elif self._dim == 4:
            stable_labels_plot = go.Scatter3d(x=x, y=y, z=z, **plot_args)

        return stable_labels_plot

    def _create_plotly_markers(self, label_uncertainties=False):
        """
        Creates stable and unstable marker plots for overlaying on the phase diagram.

        :return: Tuple of Plotly go.Scatter (or go.Scatter3d) objects in order: (
            stable markers, unstable markers)
        """

        def get_marker_props(coords, entries, stable=True):
            """ Method for getting marker locations, hovertext, and error bars
            from pd_plot_data"""
            x, y, z, texts, energies, uncertainties = [], [], [], [], [], []
            
            for coord, entry in zip(coords, entries):
                energy = round(self._pd.get_form_energy_per_atom(entry), 3)

                entry_id = getattr(entry, "entry_id", "no ID")
                comp = entry.composition

                if hasattr(entry, "original_entry"):
                    comp = entry.original_entry.composition

                formula = comp.reduced_formula
                clean_formula = self._htmlize_formula(formula)
                label = f"Comp: {clean_formula} <br>" f"Form_E: {energy} eV/atom <br>"
                    
                # This is where Jiadong makes changes
                
#                 if entry.name == self.new_entries[-1].name:
#                 if not entry.is_element:
#                     mod_entries = [
#                         e
#                         for e in self._pd._stable_entries
#                         if e.name != entry.name
#                     ]
#                     mod_cpd = PhaseDiagram(mod_entries)
# #                     print(entry.name)
#                     invE = mod_cpd.get_decomp_and_e_above_hull(entry, allow_negative=True)[1]
#                     label += f"Inv_E: {invE} <br>"
                
                label += f"{self.reactions_dict[entry.name].__str__()}"

#                 for ee in self.new_entries:
#                     if entry.name == ee.name:
#                         ind = self.new_entries.index(ee)
# #                         if ind < len(self.new_entries)-1:
#                         label += f"{self.reactions_list[ind].__str__()}"
# #                         elif ind == len(self.new_entries)-1:
# #                             label += f"{ee.name}"

                if not stable:
                    e_above_hull = round(self._pd.get_e_above_hull(entry), 3)
                    if e_above_hull > self.show_unstable:
                        continue
                    label += f" (+{e_above_hull} eV/atom)"
                    energies.append(e_above_hull)
                else:
                    uncertainty = 0
                    if (
                        hasattr(entry, "correction_uncertainty_per_atom")
                        and label_uncertainties
                    ):
                        uncertainty = round(entry.correction_uncertainty_per_atom, 4)
                        label += f"<br> (Error: +/- {uncertainty} eV/atom)"

                    uncertainties.append(uncertainty)
                    energies.append(energy)

                texts.append(label)

                x.append(coord[0])
                y.append(coord[1])

                if self._dim == 3:
                    z.append(energy)
                elif self._dim == 4:
                    z.append(coord[2])

            return {
                "x": x,
                "y": y,
                "z": z,
                "texts": texts,
                "energies": energies,
                "uncertainties": uncertainties,
            }

        stable_coords, stable_entries = (
            self.pd_plot_data[1].keys(),
            self.pd_plot_data[1].values(),
        )
        unstable_entries, unstable_coords = (
            self.pd_plot_data[2].keys(),
            self.pd_plot_data[2].values(),
        )
        
        stable_props = get_marker_props(stable_coords, stable_entries)

        unstable_props = get_marker_props(
            unstable_coords, unstable_entries, stable=False
        )

        stable_markers, unstable_markers = dict(), dict()

        if self._dim == 2:
            stable_markers = plotly_layouts["default_binary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=list(stable_props["x"]),
                    y=list(stable_props["y"]),
                    name="Stable",
                    marker= dict(color="darkgreen", size=11, line=dict(color="black", width=2)),
                    opacity=0.9,
                    hovertext=stable_props["texts"],
                    error_y=dict(
                        array=list(stable_props["uncertainties"]),
                        type="data",
                        color="gray",
                        thickness=2.5,
                        width=5,
                    ),
                )
            )

            unstable_markers = plotly_layouts["default_binary_marker_settings"].copy()
            unstable_markers.update(
                dict(
                    x=list(unstable_props["x"]),
                    y=list(unstable_props["y"]),
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=6,
                        symbol="diamond",
                    ),
                    hovertext=unstable_props["texts"],
                )
            )

        elif self._dim == 3:
            stable_markers = plotly_layouts["default_ternary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=list(stable_props["y"]),
                    y=list(stable_props["x"]),
                    z=list(stable_props["z"]),
                    name="Stable",
                    marker=dict(
                        color="black",
                        size=12,
                        opacity=0.8,
                        line=dict(color="black", width=3),
                    ),
                    hovertext=stable_props["texts"],
                    error_z=dict(
                        array=list(stable_props["uncertainties"]),
                        type="data",
                        color="darkgray",
                        width=10,
                        thickness=5,
                    ),
                )
            )

            unstable_markers = plotly_layouts["default_ternary_marker_settings"].copy()
            unstable_markers.update(
                dict(
                    x=unstable_props["y"],
                    y=unstable_props["x"],
                    z=unstable_props["z"],
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=6,
                        symbol="diamond",
                        colorbar=dict(
                            title="Energy Above Hull<br>(eV/atom)", x=0.05, len=0.75
                        ),
                    ),
                    hovertext=unstable_props["texts"],
                )
            )

        elif self._dim == 4:
            stable_markers = plotly_layouts["default_quaternary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=stable_props["x"],
                    y=stable_props["y"],
                    z=stable_props["z"],
                    name="Stable",
                    marker=dict(
                        color=stable_props["energies"],
                        colorscale=plotly_layouts["stable_markers_colorscale"],
                        size=8,
                        opacity=0.9,
                    ),
                    hovertext=stable_props["texts"],
                )
            )

            unstable_markers = plotly_layouts[
                "default_quaternary_marker_settings"
            ].copy()
            unstable_markers.update(
                dict(
                    x=unstable_props["x"],
                    y=unstable_props["y"],
                    z=unstable_props["z"],
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=5,
                        symbol="diamond",
                        colorbar=dict(
                            title="Energy Above Hull<br>(eV/atom)", x=0.05, len=0.75
                        ),
                    ),
                    hovertext=unstable_props["texts"],
                    visible="legendonly",
                )
            )

        stable_marker_plot = (
            go.Scatter(**stable_markers)
            if self._dim == 2
            else go.Scatter3d(**stable_markers)
        )
        unstable_marker_plot = (
            go.Scatter(**unstable_markers)
            if self._dim == 2
            else go.Scatter3d(**unstable_markers)
        )

        return stable_marker_plot, unstable_marker_plot

    def show(self, *args, **kwargs):
        r"""
        Draw the phase diagram using Plotly (or Matplotlib) and show it.

        Args:
            *args: Passed to get_plot.
            **kwargs: Passed to get_plot.
        """
        
        filename = kwargs.pop('filename', 'path_to_file.html')
        print("filename",filename)
        fig = self.get_plot(*args, **kwargs)
        if self.target_entry is not None:
            for coords, entry in self.pd_plot_data[1].items():
                print(entry.name)
                if entry.name == self.target_entry.name:
                    x_coord = coords[0]
                    y_coord = coords[1]
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=[x_coord],
                    y=[y_coord],
                    marker=dict(
                        color='rgba(135, 206, 250, 0.5)',
                        size=30,
                        line=dict(
                            color='MediumPurple',
                            width=4
                        )
                    ),
                    showlegend=False
                )
            )

        fig.write_html(filename)
        fig.show()
         
        
        
        
        
        
        
        
        
        
        
        
        
        