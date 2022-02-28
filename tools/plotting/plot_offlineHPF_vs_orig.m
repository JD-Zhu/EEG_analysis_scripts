figure; hold on;
plot(freq.freq, mean(freq.powspctrm));
plot(freq.freq, mean(freq_HPF.powspctrm));
xlim([0 0.2]);
xlabel('Frequency (Hz)');
ylabel('Absolute power (uV^2)');
legend({'original', 'offline 0.01Hz HPF'})
hold off;