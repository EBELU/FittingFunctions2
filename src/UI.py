from .gaussian_fitting import fit_gaussian, calibrate
from matplotlib.widgets import SpanSelector, Button, RadioButtons, TextBox
import matplotlib.pyplot as plt
import numpy as np



class UI_wrapper:
    x_data = None
    y_data = None
    UI = None
    calibrated = False
    
    def set_data(x_data, y_data):
        UI_wrapper.x_data = x_data
        UI_wrapper.y_data = y_data
    
    def run():
        UI_wrapper.UI = ff2_UI(UI_wrapper.x_data, UI_wrapper.y_data, UI_wrapper.calibrated)

class ff2_UI:
    def __init__(self, x_data, y_data, calibrated):
        self.fig, (self.ax, self.bx) = plt.subplots(1, 2, figsize = (20, 10))
        self.fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace = 0.02)
        self.peaks = []
        self.span_selector = SpanSelector(self.ax, self.on_select, direction='horizontal', 
                                          useblit=True, props=dict(alpha=0.5, facecolor='red'))
        self.x_data, self.y_data = x_data, y_data

        self.ax.plot(x_data, y_data)
        self.bx.axis("off")
        
        button_clear_axis = self.fig.add_axes([0.5, 0.85, 0.05, 0.05])


        self.extra_axes = {"b_clear" : button_clear_axis,
                           "b_print_peaks": self.fig.add_axes([0.55, 0.85, 0.05, 0.05]),
                            "b_l_corr" : self.fig.add_axes([0.6, 0.7, 0.05, 0.05]),
                            "b_r_corr" : self.fig.add_axes([0.65, 0.7, 0.05, 0.05])}
        
        self.buttons = {"clear": Button(button_clear_axis, 'Clear'),
                        "print_peaks": Button(self.extra_axes["b_print_peaks"], "Print Peaks"),
                        "l_corr" : TextBox(self.extra_axes["b_l_corr"], "Scatter \nCorr ", initial="L"),
                        "r_corr" : TextBox(self.extra_axes["b_r_corr"], "", initial="R")}
        
        self.buttons["l_corr"].on_submit(self.change_corr_wrapper("L"))
        self.buttons["r_corr"].on_submit(self.change_corr_wrapper("R"))
                
        self.peak_axes = dict([])
        self.peak_buttons = dict([])
        
        self.bx_t = self.bx.text(.15,.4, "", fontsize = 11)
        self.logger_text = self.bx.text(0,-0.1, "Log", fontsize = 7)
        
        self.buttons["clear"].on_clicked(self.clear)
        self.buttons["print_peaks"].on_clicked(self.print_peaks)
        
        if not calibrated:
            self.extra_axes.update({"b_calib" : self.fig.add_axes([0.5, 0.80, 0.05, 0.05]),
                                       "tb_calib" : self.fig.add_axes([0.60, 0.80, 0.20, 0.05])})
            
            self.tb_calib = TextBox(self.extra_axes["tb_calib"], "Energies \n [keV]", initial="E1, E2, ...")
            self.tb_calib.on_submit(self.set_calibration_data)
                        
            self.buttons.update({"calib" : Button(self.extra_axes["b_calib"], "Calibrate")})
            self.buttons["calib"].on_clicked(self.calibrate)
        
        self.log = []
        self.logger("UI initiated")

    def logger(self, text, error = False):
        if error:
            text = "ERROR!  " + text
        print(text)
        self.log.append(f"{len(self.log)}. " + str(text) + "\n")
        self.logger_text.set_text("Log:\n" + "".join(self.log[-4:]))
        plt.draw()
                
    def clear(self, arg):
        self.peaks.clear()
        self.ax.clear()
        self.ax.plot(self.x_data, self.y_data)
        self.bx_t.set_text("")
        for btn in self.peak_axes.values():
            btn.remove()
        self.peak_axes.clear()
        self.peak_buttons.clear()
        try:
            plt.draw()
        except(AttributeError):
            self.logger("AttributeError in plt.draw in clear")
        self.logger("Cleared")

    def print_peaks(self, event):
        peak_strs = []
        for peak in self.peaks:
            peak_strs.append(f"ff.fit_gaussian(x_data, y_data, region_start={round(peak._region_limits[0])}," + 
                             f"region_stop={round(peak._region_limits[1])}, corr_left={round(peak._corr_points[1])}, corr_right={round(peak._corr_points[0])}), \n")
        print(f'\u001b[42m peaks = [{"".join(peak_strs)}] \u001b[0m')

    def set_peak_button(self):
        i = len(self.peaks) - 1
        key = str(round(self.peaks[-1].mu))
        self.peak_axes.update({key: self.fig.add_axes([0.5, 0.7 - 0.05*i, 0.05, 0.05])})
        self.peak_buttons.update({key : (Button(self.peak_axes[key], key))})   
        self.peak_buttons[key].on_clicked(self.onclick_wrapper(key))
        self.logger(f"Button added: {key}")
        
    def set_radio(self):
        i = len(self.peaks) - 1
        key = str(round(self.peaks[-1].mu))
        self.bx.clear()
        self.radio = RadioButtons(self.bx, [str(round(p.mu)) for p in self.peaks])
        
    def onclick_wrapper(self, name):
        def onclick(event):
            for peak in self.peaks:
                if int(name) == round(peak.mu):
                    break
            self.bx_t.set_text(f"Region limits: {round(peak._region_limits[0])} >==< {round(peak._region_limits[1])}" +
                          "\n\n" + str(peak))
            
            self.selected = peak
            
            self.buttons["l_corr"].set_val(round(peak._corr_points[1]))
            self.buttons["r_corr"].set_val(round(peak._corr_points[0]))
            plt.draw()
        return onclick
        
    def set_calibration_data(self, event):
        Es = event.split(",")
        try:
            Es = [float(E) for E in Es]
            self.calibration_energies = Es
            self.logger("Calibration energies set")
        except:
            self.logger("Invalid energy", True)
            
    def calibrate(self, event):
        try:
            Es = self.calibration_energies
        except(AttributeError):
            self.logger("Missing calibration energies", True)
            return 0
            
        if len(self.peaks) != len(Es):
            self.logger("Number of calibration energies does not match nr of peaks", True)
            return 0
        peaks_sorted = sorted(self.peaks, key = lambda p: p.mu)
        Es_sorted = sorted(Es)
        self.logger("Peaks and energies sorted")
        X, (k,m) = calibrate(self.y_data, peaks_sorted, Es_sorted)
        
        intervals = "".join([f"({round(p[0])}, {round(p[1])})," for p in peaks_sorted])
        print(f"\u001b[42mx_data, (k,m) = ff.calibrate(y_data, [{intervals}], {Es_sorted}) \u001b[0m")
        self.logger("Calibration successful! Restarting UI")
        UI_wrapper.calibrated = True
        UI_wrapper.set_data(X, self.y_data)
        UI_wrapper.run()
        
    def change_corr_wrapper(self, side):
        def change_corr(event):
            try:
                point = float(event)
            except:
                self.logger("Invalid correction point", True)
                return
            for i, peak in enumerate(self.peaks):
                if round(peak.mu) == round(self.selected.mu):
                    break
            if side.upper() == "L":
                if round(peak._corr_points[1]) == round(point):
                    return
                new_peak = fit_gaussian(self.x_data, self.y_data, self.selected._region_limits[0], self.selected._region_limits[1],
                                       corr_left=point, corr_right=self.selected._corr_points[0])
            else:
                if round(peak._corr_points[0]) == round(point):
                    return
                new_peak = fit_gaussian(self.x_data, self.y_data, self.selected._region_limits[0], self.selected._region_limits[1],
                                           corr_right=point, corr_left=self.selected._corr_points[1])
            self.peaks[i] = new_peak
            self.logger(f"Peak {round(peak.mu)}: updated scatter corr")
            self.ax.clear()
            self.ax.plot(self.x_data, self.y_data)
            self.plot_peaks()
        return change_corr
        
    def on_select(self, region_start, region_stop):
        try:
            peak = fit_gaussian(self.x_data, self.y_data, region_start, region_stop)
            if peak:
                self.logger(f"Peak fitted: {round(peak.mu)}")
                self.peaks.append(peak)
                self.plot_peaks()
                self.set_peak_button()
        except:
            self.logger("Fitting failed!", True)
            pass

    
    def plot_peaks(self):
        for peak in self.peaks:
            self.bx_t.set_text(f"Region limits: {round(peak._region_limits[0])} >==< {round(peak._region_limits[1])}" +
                          "\n\n" + str(peak))
            self.ax.plot(peak._x, peak._y, label="Data", color="deepskyblue")
            self.ax.plot(peak._x, peak._corr_f(peak._x) + peak.value(peak._x),
                     label="Fitted Gaussian", color="indigo")
            if peak._corr_points != [None, None]:
                try:
                    self.ax.plot(peak._x, peak._corr_f(peak._x), ":",
                             label="Scatter correction", color="red")
                    self.ax.vlines(peak._corr_points[0], 0, (peak._corr_f(
                        peak.mu)+peak.A) * 1, color="red")
                    self.ax.vlines(peak._corr_points[1], 0, (peak._corr_f(
                        peak.mu)+peak.A) * 1, color="red")
                except(AttributeError):
                    pass

        try:
            plt.draw()
        except(AttributeError):
            self.logger("AttributeError in plt.draw in plot_peaks")
       
