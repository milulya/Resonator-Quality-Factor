import numpy as np
from os import path
from scipy import linalg, optimize
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d
from generate import DataGenerator as DataGen
plt.style.use('ggplot')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


class Resonator:
    def __init__(self):
        self.measurements_path = {}
        self.measurements = {}
        self.measurements_names = np.array([])
        self.Qls = {}
        self.Qcs = {}
        self.Qis = {}
        self.asymmetrys = {}
        self.number_of_measurements = 0

    def new_measurement(self, msrmnt_path, name):
        self.measurements_path[name] = msrmnt_path
        self.measurements[name] = Measurement(msrmnt_path)
        self.measurements[name].name = name
        self.measurements[name].measure(plot_data=False, verbose=False)
        self.measurements[name].measurememt_number = self.number_of_measurements
        self.number_of_measurements = self.number_of_measurements + 1
        self.Qls[name] = self.measurements[name].Ql
        self.Qcs[name] = np.abs(self.measurements[name].Qc)
        self.Qis[name] = self.measurements[name].Qi
        self.asymmetrys[name] = self.measurements[name].asymmetry
        self.measurements_names = np.append(self.measurements_names, name)

    def plot_values_dict(self, val_dict, ax=None, plot_type='plot', yscale='linear', **kwargs):
        return_vals = False
        if ax is None:
            return_vals = True
            fig, ax = plt.subplots()

        try:
            # converitng measuremmnts names to float
            x_values = self.measurements_names.astype('float')
            sorted_x_args = x_values.argsort()
            y_values = [val_dict[name_] for name_ in self.measurements_names[sorted_x_args]]
        except ValueError:  # couldn't convert names to float
            # using measurements number as x values
            x_values = []
            y_values = []
            for measurement in self.measurements:
                x_values.append(measurement.measurememt_number)
                y_values.append(val_dict[measurement.name])

            sorted_x_args = x_values.argsort()
            y_values = y_values[sorted_x_args]
            x_values = x_values[sorted_x_args]

        if plot_type == 'plot':
            ax.plot(x_values[sorted_x_args], y_values, **kwargs)
        elif plot_type == 'scatter':
            ax.scatter(x_values[sorted_x_args], y_values, **kwargs)
        else:
            raise ValueError('plot_type argument is invalid')

        ax.set_yscale(yscale)

        if return_vals is True:
            return fig, ax

    # def plot_measurement_parameter(self):

    def get_msrmnt_atrr(self, attribute):

        attributes = {} #[[name_, None] for name_ in self.measurements_names]
        for i, name_ in enumerate(self.measurements_names):
            attributes[name_] = getattr(self.measurements[name_], attribute)
        return attributes






class Measurement:

    def __init__(self, vna_csv_data_path=None, data_tuple=None, config = 'T', s_mat_element = '21'):
        if not (vna_csv_data_path is None):
            self.data_path = path.normpath(vna_csv_data_path)
            self.data_dictionary = self.vnacsvreader()
            self.data_matrix_raw = self.dict2mat(s_mat_element=s_mat_element, data_dict=self.data_dictionary)
            self.frequencies_raw = self.data_matrix_raw[:, 0]
            self.frequencies = self.frequencies_raw
            self.z_data_raw = self.data_matrix_raw[:, 1] + 1j*self.data_matrix_raw[:, 2]
        elif not (data_tuple is None):
            self.frequencies_raw = data_tuple[0]
            self.frequencies = self.frequencies_raw
            self.z_data_raw = data_tuple[1]
        else:
            raise(ValueError('no data has given, supply path to csv or tuple in the format: (frequencies, z_data)'))
        if config == 'T' or config.lower()=="circulator":
            self.config = config
        else:
            rasie(ValueError('setup configuration is not clear, specift "T" or "circulator"'))
        self.z_data_undelayed = None
        self.z_data_calibrated = None
        self.z_data_generated = None
        self.circle = None
        self.calibrated_circle = None
        self.delay_rough_estimation = None
        self.Ql = None
        self.Qc = None
        self.Qi = None
        self.Qi_no_correction = None
        self.delay = None
        self.env_factor = None
        self.fr = None
        self.theta_0 = None
        self.asymmetry = None
        self.kappa = None
        self.mse = None
        self.name = None

    def measure(self, plot_data=False, verbose=False):
        # self.plot_measured_data()
        self.delay_rough_estimation, poly = self.estimate_cable_delay()
        self.delay = self.calculate_cable_delay(self.frequencies, self.z_data_raw)
        self.z_data_undelayed = self.correctdelay(self.frequencies, self.z_data_raw, self.delay)
        # fitting circle to data after delay correction
        xc, yc, r0 = self.fit_circle(self.z_data_undelayed)
        self.circle = {'xc': xc, 'yc': yc, 'r0': r0}

        self.fr, self.theta_0, self.Ql, origin_phase = self.fit_phase(self.frequencies, self.z_data_undelayed)
        self.env_factor = self.calc_env_factor()
        self.z_data_calibrated, self.calibrated_circle = self.calibrate()
        self.Qc, self.Qi = self.calculate_quality_factors()
        # generating synthetic data according to the model theory using parameters calculated in the above code
        gen = DataGen(self.Ql, np.abs(self.Qc), self.asymmetry, self.frequencies, self.fr, np.abs(self.env_factor),
                      np.angle(self.env_factor), self.delay, config=self.config)
        self.z_data_generated = gen.z_data_env

        # calculating erros
        self.calculate_errors()

        if plot_data is True:
            # phase fit plot
            plt.figure()
            plt.plot(self.frequencies, origin_phase, label='Phase')
            plt.plot(self.frequencies, self.theta_0 + 2. * np.arctan(2. * self.Ql * (1 - self.frequencies / self.fr)),
                     label='fit')
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Phase [Rad]')
            plt.suptitle('Fit of phase response\nfor the resonance circle translated to the origin')
            plt.legend()
            plt.axvline(self.fr, color='black', alpha=0.5, ls='--', label='resonance')
            # cable delay correction plot
            plt.figure()
            plt.scatter(self.z_data_undelayed.real, self.z_data_undelayed.imag, label='Delay Corrected Data')
            theta = np.linspace(0, 2 * np.pi, len(self.frequencies))
            fitted_circle = r0 * np.exp(1j * theta) + np.complex(xc, yc)
            plt.scatter(fitted_circle.real, fitted_circle.imag, s=8, label='Circle Fit')
            plt.xlabel('Real')
            plt.ylabel('Imag')
            plt.suptitle(f'circle fit after calculating delay using list squares\ncable delay:{self.delay:.3E}')
            plt.legend()
            # calibration plot
            self.plot_calibration()
            # data generated from calculated values plotted over measured data
            self.sanity_plot()

        if verbose is True:
            print(f'Total Quality factor is {self.Ql:.3E}\n'
                  f'Absolute Value of coupling quality factor is {np.abs(self.Qc):.3E}\n'
                  f'Real Part of coupling Quality factor Qc is {self.Qc.real:.3E}\n'
                  f'Intrinsic Quality Factor is {self.Qi:.3E}\n'
                  f'Intrinsic Quality Factor without assymtery correction is {self.Qi_no_correction:.3E}\n'
                  f'Resonance Frequency is {self.fr:.5E}')

        return self.Ql, self.Qc, self.Qi, self.fr

    def calculate_errors(self):
        # calculating Mean Squared Error of the generated data compared to the measured data

        # normalizing the data
        xc, yc, r0 = self.fit_circle(self.z_data_raw)
        z_data_raw = (self.z_data_raw - np.complex(xc, yc))/r0
        z_data_generated = (self.z_data_generated - np.complex(xc, yc))/r0
        mse_x = np.sum((z_data_generated.real - z_data_raw.real)**2) / len(self.z_data_raw)
        mse_y = np.sum(z_data_generated.imag - z_data_raw.imag)**2 / len(self.z_data_raw)
        self.mse = mse_x + mse_y
        return self.mse

        # TODO: Calculate SNR

    def calculate_quality_factors(self):
        if self.config == 'T':
            abs_Qc = self.Ql / (2*self.calibrated_circle['r0'])
        elif self.config == 'circulator':
            abs_Qc = self.Ql / (self.calibrated_circle['r0'])
        # calculating the asymmetry phi
        # shifting the off resonance point to the origin
        xc = self.calibrated_circle['xc'] - 1
        yc = self.calibrated_circle['yc']
        # angle_Qc = -np.arcsin(self.calibrated_circle['yc'] / self.calibrated_circle['r0'])

        if xc > 0:
            if yc > 0:
                angle_Qc = -np.pi + np.arcsin(yc / self.calibrated_circle['r0'])
            elif yc <= 0:
                angle_Qc = np.pi + np.arcsin(yc / self.calibrated_circle['r0'])
        elif xc <= 0:
            angle_Qc = -np.arcsin(self.calibrated_circle['yc'] / self.calibrated_circle['r0'])

        Qc = abs_Qc*np.exp(-1j*angle_Qc)
        self.asymmetry = angle_Qc
        Qi = 1 / (1/self.Ql - (1/Qc).real)
        self.Qi_no_correction = 1 / (1/self.Ql - 1/np.abs(Qc))
        return Qc, Qi

    def calibrate(self):
        z_data_calibrated = self.z_data_undelayed/self.env_factor
        circle_center = np.complex(self.circle['xc'], self.circle['yc'])
        calibrated_circle_center = circle_center/self.env_factor
        calibrated_circle_radius = self.circle['r0']/np.abs(self.env_factor)
        calibrated_circle = {'xc': calibrated_circle_center.real,
                             'yc': calibrated_circle_center.imag,
                             'r0': calibrated_circle_radius}
        return z_data_calibrated, calibrated_circle

    def plot_measured_data(self):
        fig_mag_phase, ax_mag_phase = plt.subplots(2, 1)
        frequencies = self.frequencies
        z_mag = np.abs(self.z_data_raw)
        z_phase = np.angle(self.z_data_raw)

        ax_mag_phase[0].plot(frequencies, z_mag, label='Measured')
        ax_mag_phase[0].set_ylabel('Magnitude')
        ax_mag_phase[0].set_title('S21 Magnitude Vs Frequency')
        ax_mag_phase[1].plot(frequencies, z_phase, label='Measured')
        ax_mag_phase[1].set_title('S21 Phase Vs Frequency')
        ax_mag_phase[1].set_ylabel('Phase [Rad]')
        for ax in ax_mag_phase:
            ax.set_xlabel('Frequency [Hz]')
        fig_circle, ax_circle = plt.subplots()
        self.plot_circle(ax_circle, self.z_data_raw, label='Measured')
        ax_circle.set_title('Polar Plot')
        fig_circle.suptitle('Measured Data')
        fig_mag_phase.suptitle('Measured Data')

        return fig_mag_phase, fig_circle, ax_mag_phase, ax_circle

    def sanity_plot(self):
        """
        this function plots genreated data using the calculated values on top of the measured data
        :return:
        """
        fig_mag_phase, fig_circle, ax_mag_phase, ax_circle = self.plot_measured_data()
        ax_mag_phase[0].plot(self.frequencies, np.abs(self.z_data_generated), label='Generated')
        ax_mag_phase[0].axvline(self.fr, color='black', alpha=0.5, ls='--')
        ax_mag_phase[0].legend()
        ax_mag_phase[1].plot(self.frequencies, (np.angle(self.z_data_generated)), label='Generated')
        ax_mag_phase[1].axvline(self.fr, color='black', alpha=0.5, ls='--')
        ax_mag_phase[1].legend()
        self.plot_circle(ax_circle, self.z_data_generated, s=8, label='Generated')
        ax_circle.legend()
        fig_mag_phase.suptitle('sanity plots')
        fig_circle.suptitle('sanity plot')

    @staticmethod
    def plot_circle(ax, z_data, **kwargs):
        """
        plots a circle on and axis. taking care of the axis limits.
        the function can plots many circles on the same axis while taking care of the axis limits and lines
        :param ax: axis to plot
        :param z_data: circle complex data
        :param kwargs: kew word arguments to pass to the scatter function
        :return:
        """

        ax.scatter(z_data.real, z_data.imag, **kwargs)
        x_lim = np.max(ax.get_xlim())
        y_lim = np.max(ax.get_ylim())
        data_xlim = 1.1*np.max(np.abs(z_data.real))
        data_ylim = 1.1*np.max(np.abs(z_data.imag))
        data_xylim = np.max([data_xlim, data_ylim])
        if data_xylim > x_lim:
            ax.set_xlim([-data_xylim, data_xylim])
        if data_xylim > y_lim:
            ax.set_ylim([-data_xylim, data_xylim])

        lines = ax.get_lines()
        if len(lines) == 0:
            plt.axvline(0, color='black', alpha=0.5, ls='--')
            plt.axhline(0, color='black', alpha=0.5, ls='--')
        ax.set_xlabel('Real')
        ax.set_ylabel('Imaginary')

    def plot_calibration(self):
        xc = self.circle['xc']
        yc = self.circle['yc']
        r0 = self.circle['r0']
        plt.figure()
        ax = plt.gca()
        self.plot_circle(ax, self.z_data_undelayed, label='Delay Corrected Data')
        theta = np.linspace(0, 2 * np.pi, len(self.frequencies))
        circle = r0 * np.exp(1j * theta) + np.complex(xc, yc)
        res_point = r0*np.exp(1j*self.theta_0) + np.complex(xc, yc)
        beta = self.periodic_angle(self.theta_0 + np.pi)
        off_res_point = r0*np.exp(1j*beta) + np.complex(xc, yc)
        self.plot_circle(ax, circle, s=8, label='Circle Fit')
        env_factor = self.env_factor
        calibrated_circle = circle / env_factor
        self.plot_circle(ax, calibrated_circle, s=8, c=colors[4], label=' Calibrated Circle Fit')
        plt.scatter(res_point.real, res_point.imag, marker='*', s=64, c='fuchsia', label='resonance point')
        plt.scatter(off_res_point.real, off_res_point.imag, marker='*', s=64, c='lime', label='off-resonance point')
        plt.scatter((res_point/env_factor).real, (res_point/env_factor).imag, s=64, c='fuchsia', marker='*')
        plt.scatter((off_res_point/env_factor).real, (off_res_point/env_factor).imag, s=64, c='lime', marker='*')
        plt.legend()

    def calc_env_factor(self):
        """
        calculating the environment functor multiplying the S21 canonical form
        :return:
        """
        if (self.circle is None) or (self.z_data_undelayed is None):
            raise ValueError('fitted circle and delay corrected data are missing,'
                             ' Run "fit_phase" before running "calc_env_factor')
        xc = self.circle['xc']
        yc = self.circle['yc']
        r0 = self.circle['r0']
        beta = self.periodic_angle(self.theta_0 + np.pi)
        off_resonance_point = np.complex(xc + r0 * np.cos(beta), yc + r0 * np.sin(beta))
        # since off resonance point should be 1, the environment factor is equal to the off resonanc epoint value
        self.env_factor = off_resonance_point
        return off_resonance_point

    def fit_phase(self, frequencies=None, z_data=None, plot=False):
        """
        this function fits the phase response of S21 circle shifted to the center
        extracting accurate values for  resonance frequency fr, total Quality factor QL and
        of resonance point angle theta_0.

        :param frequencies:
        :param z_data:
        :param plot:
        :return:
        """
        if (frequencies is None) and (z_data is None):
            frequencies = self.frequencies
            z_data = self.z_data_undelayed
        elif (frequencies is None) ^ (z_data is None):
            raise ValueError('please supply both frequencies and z_data or none of them')

        if self.delay is None:
            raise NameError('missing cable delay value')

        # z_data = self.correctdelay(frequencies, z_data, self.delay)
        # self.z_data_undelayed = z_data.copy()
        xc, yc, r0 = self.fit_circle(z_data)
        # self.circle = {'xc': xc, 'yc': yc, 'r0': r0}
        # xc = self.circle['xc']
        # yc = self.circle['yc']
        z_data = z_data - np.complex(xc, yc)
        phase = np.unwrap(np.angle(z_data))

        # calculating initial value for resonance frequency
        phase_smooth = gaussian_filter1d(phase, 30)
        phase_derivative = np.gradient(phase_smooth)
        fr_initial_val = frequencies[np.argmax(phase_derivative)]

        # calculating initial value for Ql
        try:
            freq_truncated, z_truncated = self.narrow_band(width_multiplier=1)
            delta_f = freq_truncated[-1] - freq_truncated[0]
        except ValueError:
            delta_f = frequencies[-1] - frequencies[0]

        Ql_intial_value = fr_initial_val / delta_f

        theta_intial_value = phase[0] - np.pi
        initial_values = [fr_initial_val, theta_intial_value, Ql_intial_value]

        # defining the residuals function  for least squares fit
        def residuals(fit_parametrs, frequencies, phase):
            theta = fit_parametrs[1]
            Ql = fit_parametrs[2]
            fr = fit_parametrs[0]
            res = theta + 2.*np.arctan(2.*Ql*(1 - frequencies/fr)) - phase
            return res

        # making first a fit for each parameter separately, for better convergence
        for i, parameter in enumerate(initial_values):
            def residuals_partial(parameter, frequencies, phase):
                initial_values[i] = parameter
                return residuals(initial_values, frequencies, phase)

            parameter_semi_optimized = optimize.least_squares(residuals_partial, parameter, args=(frequencies, phase),
                                                              ftol=1e-12, xtol=1e-12)
            initial_values[i] = parameter_semi_optimized.x[0]
        # final optimization for all parameters together starting with the values of independently fitted parameters
        optimized = optimize.least_squares(residuals, initial_values, args=(frequencies, phase), ftol=1e-12, xtol=1e-12)
        calculated_values = optimized['x']
        fr, theta_0, Ql = calculated_values
        self.fr = fr
        self.theta_0 = theta_0
        self.Ql = Ql
        return fr, theta_0, Ql, phase

    @staticmethod
    def fit_circle(z_data, plot=False):
        """
        this function fits circle according to the algorithm described in:

        Probst, S.;. F. B.;. P. A.;. A. V.;. M.
        “Efficient and Robust Analysis of Complex Scattering Data Under Noise in Microwave Resonators.”
        Review of Scientific -Instruments, vol. 86, no. 2, AIP Publishing LLC, Feb. 2015, p. 024706,
          doi:10.1063/1.4907935.

        :return:
        """
        # Normalizing data point for higher numerical accuracy
        x_avg = (np.max(z_data.real) + np.min(z_data.real)) / 2
        y_avg = (np.max(z_data.imag) + np.min(z_data.imag)) / 2
        z_data_shifted = z_data - (x_avg + 1j * y_avg)
        z_norm = np.max(np.abs(z_data_shifted))
        z_data_normalized = z_data_shifted / z_norm

        x = z_data_normalized.real
        y = z_data_normalized.imag
        z = np.power(x, 2) + np.power(y, 2)

        # calculating the matrix of moments
        moments = np.zeros((4, 4))
        for row, array1 in enumerate([z, x, y, np.ones(len(z_data))]):
            for column, array2 in enumerate([z, x, y, np.ones(len(z_data))]):
                moments[row, column] = np.sum(array1*array2)
        constrain_matrix = np.array([[0, 0, 0, -2],
                                     [0, 1, 0, 0],
                                     [0, 0, 1, 0],
                                     [-2, 0, 0, 0]])
        [eigen_vals, eigen_vecs] = linalg.eig(moments, constrain_matrix, right=True)
        # verifying there are only real eigen values
        if not np.all(np.isreal(eigen_vals)):
            raise ValueError('imaginary eigen values calculated, validate data')
        else:
            eigen_vals = eigen_vals.real
        # finding the smallest non-negative eigen value and th corresponding eigen vector
        eigen_vecs = eigen_vecs[:, eigen_vals >= 0]
        eigen_vals = eigen_vals[eigen_vals >= 0]
        eta = eigen_vals[eigen_vals.argmin()]
        A_vec = eigen_vecs[:, eigen_vals.argmin()]

        xc = -A_vec[1] / (2. * A_vec[0])
        yc = -A_vec[2] / (2. * A_vec[0])
        # The term *sqrt term corrects for the constraint, because it may be
        # altered due to numerical inaccuracies during calculation
        r0 = (1. / (2. * np.absolute(A_vec[0]))) * np.sqrt(
            A_vec[1] * A_vec[1] + A_vec[2] * A_vec[2] - 4. * A_vec[0] * A_vec[3])

        if plot:
            theta = np.linspace(0, 2*np.pi, 1000)
            circle = r0*z_norm*np.exp(1j*theta) + xc*z_norm + x_avg +1j*(yc*z_norm + y_avg)
            plt.figure()
            plt.scatter(z_data.real, z_data.imag, label='data')
            plt.scatter(circle.real, circle.imag, label='fit')
            plt.legend()
        return xc*z_norm + x_avg, yc*z_norm + y_avg, r0*z_norm

    def calculate_cable_delay(self, frequencies=None, z_data=None):
        """
        calculating cable delay using least square fit of the data for a circle for different values of delay
        :param z_data:
        :return:
        """
        if z_data is None:
            z_data = self.z_data_raw
            frequencies = self.frequencies
        # first part - clculate cable delay by minimizing the deviance from a shpae of a circle
        def residuals(delay):
            z_data_ = self.correctdelay(frequencies, z_data, delay[0])
            xc, yc, r0 = self.fit_circle(z_data_)
            res = np.sqrt((z_data_.real - xc)**2 + (z_data_.imag - yc)**2) - r0
            # normalizing residulas for beeter convergence
            res = res / r0
            return res
        delay_upper_bound = 100e-9
        if self.delay_rough_estimation <= delay_upper_bound:
            initial_guess = self.delay_rough_estimation
        else:
            initial_guess = delay_upper_bound

        optimized = optimize.least_squares(residuals, initial_guess/2, bounds=(0, delay_upper_bound), xtol=1e-12,
                                           ftol=1e-12, gtol=1e-12)
        cable_delay = optimized.x[0]
        # second part - fine adjumnets using phase response curve
        z_data_undelayed = self.correctdelay(frequencies, z_data, cable_delay)
        fr, theta_0, Ql, origin_phase = self.fit_phase(frequencies, z_data_undelayed)



        self.delay = cable_delay
        return cable_delay

    @staticmethod
    def correctdelay(frequencies, z_data, delay):
        delay_correction = np.exp(1j * 2. * np.pi * frequencies * delay)
        z_data_corrected = z_data * delay_correction
        return z_data_corrected

    def estimate_cable_delay(self, frequencies=None, z_data=None, plot=False):
        """
        fits the phase of the data to linear function in order to give rough estimate for the cable delay
        :param data_matrix: first column is frequency second colums is real and third is imaginary
        :return: cable delay in nano second
        """
        if (frequencies is None) and (z_data is None):
            frequencies = self.frequencies
            z_data = self.z_data_raw

        phase = np.unwrap(np.angle(z_data))
        poly = np.polyfit(frequencies,phase, deg=1)
        cable_delay = np.abs(poly[0]/(2.*np.pi))
        if self.delay_rough_estimation is None:
            self.delay_rough_estimation = cable_delay

        return cable_delay, poly

    @staticmethod
    def periodic_angle(angle):
        """
        maps an angle to the interval [-pi pi)

        """
        return (angle + np.pi) % (2*np.pi) - np.pi

    def dict2mat(self, s_mat_element='21', data_dict=None):
        if data_dict is None:
            data_dict = self.data_dictionary
        keys = data_dict.keys()
        freqs = data_dict['Freq(Hz)']
        for key in keys:
            if ('REAL' in key) and (s_mat_element in key):
                real = data_dict[key]
            elif ('IMAG' in key) and (s_mat_element in key):
                imag = data_dict[key]
            else:
                continue
        try:
            data_mat = np.column_stack((freqs, real, imag))
        except UnboundLocalError as err:
            print('data file is not matched with expected format (REAL IMAG or s matrix element')
            raise err

        return data_mat

    def vnacsvreader(self, csvpath=None, skiprows=6, droprows=2):
        """
        this function parse csv files created by the keysight network analyzer
        :param csvpath: path to csv file
        :param skiprows: number of rows to skip at start (not including the columns titles)
        :param droprows: rows to drop at the end of the file
        :return: a dictionary with the values as specified in the columns titles
        dictionary values:
            """
        # TODO add support for both kinds of data, magnitude & phase, real & imaginary
        if csvpath is None:
            csvpath = self.data_path
        else:
            csvpath = path.normpath(csvpath)
        with open(csvpath, 'r') as fp:
            lines = fp.readlines()
            headers = lines[skiprows][:-1].split(',')
        lines = lines[skiprows + 1:-droprows]
        lines = [line[:-1].split(',') for line in lines]
        nparray = np.transpose(np.array(lines, dtype='float'))
        dictionary = dict((zip(headers, nparray)))
        return dictionary

    def narrow_band(self, width_multiplier=1):
        """
        this function truncate the data so all the oints will be closed enough to resosnace.
        it doing this by doing a rough estimation of the width |S21| line shape and takes the data to multipliers
        of this width
        :param width_multiplier: multipliers of the line width, 4 should be good
        :param plot: if TRue will generate plot of scatter x-y
        :return: truncated frequencies & z_data
        """
        frequencies = self.frequencies
        z_data = self.z_data_raw

        rough_delay, _ = self.estimate_cable_delay()
        z_data_corrected = self.correctdelay(frequencies, z_data, rough_delay)

        phase_smooth = gaussian_filter1d(np.unwrap(np.angle(z_data_corrected)), 30)
        phase_derivative = np.gradient(phase_smooth)
        phase_derivative_smooth = gaussian_filter1d(phase_derivative, 80)
        phase_second_derivative = gaussian_filter1d(np.gradient(phase_derivative_smooth), 50)

        band_frequencies = frequencies[phase_second_derivative.argmax(): phase_second_derivative.argmin()]
        if len(band_frequencies) < 2:
            raise (ValueError('was not able to truncate data'))
        bandwidth = band_frequencies[-1] - band_frequencies[0]
        extension_multiplier = (width_multiplier-1)/2
        lowest_freq = band_frequencies[0] - bandwidth*extension_multiplier
        highest_freq = band_frequencies[-1] + bandwidth*extension_multiplier

        truncated_freq = frequencies[np.logical_and(lowest_freq < frequencies, frequencies < highest_freq)]
        truncated_z_data = z_data[np.logical_and(lowest_freq < frequencies, frequencies < highest_freq)]

        return truncated_freq, truncated_z_data

    def smoth_mag_noise(self, sigma=5):
        magnitude = gaussian_filter1d(np.abs(self.z_data_raw), sigma)
        angle = gaussian_filter1d(np.angle(self.z_data_raw), sigma)
        self.z_data_raw = magnitude*np.exp(1j*angle)
