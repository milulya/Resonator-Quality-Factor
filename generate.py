import numpy as np
import matplotlib.pyplot as plt


class DataGenerator:

    def __init__(self, Ql, Qc, phi, frequencies, fr, a, alpha, delay, config='T'):
        self.Ql = Ql
        self.Qc = Qc
        self.phi = phi
        self.fr = fr
        self.frequencies = frequencies
        self.a = a
        self.alpha = alpha
        self.delay = delay
        self.z_data_circle = None
        self.z_data_env = None
        self.z_data_env_noise = None
        self.z_data_noise = None
        if config == 'T' or config == "circulator":
            self.config = config
        else:
            rasie(ValueError('setup configuration is not clear, specift "T" or "circulator"'))
        self.resonance_circle(self.config)
        self.add_noise()
        self.add_environment()


    def resonance_circle(self, config):
        if not (self.z_data_circle is None):
            raise(AttributeError('data already exist, reset first'))
        if config == 'T':
            S21 = 1 - (self.Ql / self.Qc) / (1 + 2 *1j*self.Ql * (self.frequencies / self.fr - 1)) * np.exp(
                1j*self.phi)
        elif config == 'circulator':
            S21 = 1 - (2*self.Ql / self.Qc) / (1 + 2 * 1j*self.Ql * (self.frequencies / self.fr - 1)) * np.exp(
                1j*self.phi)
        self.z_data_circle = S21
        return S21

    def add_environment(self):
        if not (self.z_data_env is None):
            raise(AttributeError('data already exist, reset first'))
        environment = self.a*np.exp(1j*self.alpha)*np.exp(-1j * 2*np.pi * self.frequencies * self.delay)
        self.z_data_env = environment*self.z_data_circle
        self.z_data_env_noise = environment*self.z_data_noise
        return environment*self.z_data_noise

    def add_noise(self):
        z_data = self.z_data_circle
        r0 = np.abs(1 - self.Ql/(2*self.Qc))
        noise_size = np.random.normal(1, 0.0*r0, len(z_data))
        noise_phase = np.random.normal(0, 0.01*np.pi, len(z_data))
        noised = z_data * noise_size * np.exp(1j*noise_phase)
        self.z_data_noise = noised
        return noised

    def reset_data(self):
        self.z_data_circle = None
        self.z_data_env = None
        self.z_data_env_noise = None
        self.z_data_noise = None


if __name__ == "__main__":
    Ql = 900
    phi = 0.1
    Qc = 1000
    fr = 4.438e9
    f_range = (4.43e9, 4.45e9)
    num_points = 6401
    frequencies = np.linspace(f_range[0], f_range[1], num_points)
    a = 0.5
    alpha = 2
    delay = 10e-9
    # frequencies = frequencies[frequencies<fr]

    data_obj = DataGenerator(Ql, Qc, phi, frequencies, fr, a, alpha, delay )
    data_obj.generate()
    S21_original = data_obj.z_data_circle
    S21_env = data_obj.z_data_env
    S21_env_noise = data_obj.z_data_env_noise
    S21_noise = data_obj.z_data_noise

    fig, ax = plt.subplots(2, 2)
    ax[0, 0].scatter(S21_original.real, S21_original.imag)
    ax[0, 1].scatter(S21_noise.real, S21_noise.imag)
    ax[1, 0].scatter(S21_env.real, S21_env.imag)
    ax[1, 1].scatter(S21_env_noise.real, S21_env_noise.imag)


    for data in [S21_original, S21_env_noise]:
        fig, ax = plt.subplots(3, 1)
        ax[0].scatter(data.real, data.imag)
        ax[0].set_title('real Vs imag')
        ax[1].plot(frequencies, np.abs(data)**2)
        ax[1].set_title('abs')
        ax[2].plot(frequencies, np.unwrap(np.angle(data)))
        ax[2].set_title('phase')

    plt.show()