function numerical_temporal_analysis()
clear all;

% D = 0.0247; % um^2/sec
% D = 0.43;
D = 0.005:0.01:0.43; % um^2/sec
% D = [0.005, 0.01, 0.05, 0.1, 0.2, 0.43];
% D = [0.0123, 0.0247, 0.0378, 0.0494, 0.108, 0.216, 0.332, 0.432];

% D = 0.005:0.04:0.43; % um^2/sec
k = 1.49 * 2*pi./0.550; %um
tc = 1./(4*k^2*D);  % seconds
dt = 1/1000; % temporal spacing in seconds
num_t = 10; % total time in seconds
t = 0:dt:num_t; % seconds
t2 = [t(end:-1:2) ,t];

% figure; hold on;

% acf = exp(-t/tc);
%  acf = [acf(end:-1:2) ,acf];
%  find_correlation_decay(t, acf)

%Test diffusion coefficents with exposure time
list3 = zeros(length(D),4);
for d = 1 : length(D)
    acf = exp(-t/tc(d));
    
    %     sigma_n = sqrt(acf(1));
    % acf    = acf - acf(end);
    % acf    = acf./acf(1).*sigma_n.^2;
    
    acf = [acf(end:-1:2) ,acf];
    exposure_time = 32/1000;
    sampling = 32;
    t_s = t(1:sampling:end);
    %Testing sampling at a single D (no exposure time)
    %     sampling = 1:1:1000;
    %     list3 = zeros(length(sampling),2);
    %     for d = 1 : length(sampling)
    %         acf = exp(-t/tc);
    
    
    %% Testing Exposure time on theory curve
    %     tri = @(t,T) (T-t) .* (t<=T);
    %     E = tri(t,exposure_time);
    %     E = [E(end:-1:2), E];
    %     Enorm = E/sum(E);
    %     acf = [acf(end:-1:2) ,acf];
    %     acf_smooth = xcorr(Enorm, acf,'none');
    %     plot(acf((length(acf)+1)/2:(length(acf)+1)/2+500));hold on;
    %     plot(acf_smooth((length(acf_smooth)+1)/2:(length(acf_smooth)+1)/2+500));
    % list3(d,:) = [...
    % (4*k^2*(find_correlation_decay_set_range(t, acf))^-1)^-1,...
    % (4*k^2*(find_correlation_decay_set_range(t, acf_smooth))^-1)^-1];
    % end
    
    
    %     figure(3);
    %     plot(t, acf((length(acf)+1)/2:end));
    %         plot(acf((length(acf)+1)/2:end));
    % hold on;
    % plot(mean(list(:,(length(out)+1)/2:end)));
    
    
    n = 10; % number of signals to average 25
    list2 = zeros(n, 2);
    
    for i = 1 : n
        
        nn = 250;
        list4 = zeros(n, floor(length(acf)/sampling));
        %         list4 = zeros(n, length(acf)-32);
        for ii = 1 : nn
            signal = genarate_random_1D_signal(acf, length(acf));
            signal = sample_and_exposure_signal(signal, 32);
            %             signal = sample_signal(signal, 32);
            % signal = exposure_signal(signal, 32);
            %             signal_nonrandom = genarate_nonrandom_1D_signal(acf, length(acf));
            %              signal = smooth(signal, exposure_time/dt);
            %             signal = signal(1:sampling:end);
            %             sampled_signal_nonrandom = signal_nonrandom(1:sampling:end);
            
            
            %%% Test Sampling on non_random signal
            %             listx = zeros(1, 32);
            %             for y = 1 : sampling
            %                 sampled_signal_nonrandom = signal_nonrandom(y:sampling:end);
            %                 F_sampled_nonrandom = fftshift(fftn(sampled_signal_nonrandom));
            %                 out_sampled_nonrandom = fftshift(ifftn(ifftshift(F_sampled_nonrandom.*conj(F_sampled_nonrandom))))/numel(F_sampled_nonrandom);
            %                 listx(y) = (4*k^2*(find_correlation_decay_set_range(t2(y:sampling:end), out_sampled_nonrandom, [2, 25]))^-1)^-1;
            %             end
            
            
            
            %%% Triangular function for convolving with ACF
            %             expTime = 32; % ms
            %                     tri = @(t,T) (T-t) .* (t<=T);
            %                     E = tri(t,exposure_time);
            % E = zeros(1, length(t));
            % E(1:(exposure_time*1000/2)) = 1;
            %                     E = [E(end:-1:2), E];
            %                     Enorm = E/sum(E);
            %                     fft_Enorm = fftshift(fft(ifftshift(Enorm)));
            %         Bexp1 = sum(out .* Enorm)
            %             Bexp = xcorr(Enorm, acf,'none');
            
            % [T,S] = deconvolution(smoothed_signal',Enorm');
            
            
            %Calculate of the ACF of the signal
            F = fftshift(fftn(signal));
            % %             F_smooth = fftshift(fftn(smoothed_signal));
            %             F_sampled = fftshift(fftn(sampled_signal));
            %             F_sampled_nonrandom = fftshift(fftn(sampled_signal_nonrandom));
            %             F_nonrandom = fftshift(fftn(signal_nonrandom));
            
            
            
            out = fftshift(ifftn(ifftshift(F.*conj(F))))/numel(F);
            %             out_smooth = fftshift(ifftn(ifftshift(F_smooth.*conj(F_smooth))))/numel(F_smooth);
            %             out_sampled = fftshift(ifftn(ifftshift(F_sampled.*conj(F_sampled))))/numel(F_sampled);
            %             out_nonrandom = fftshift(ifftn(ifftshift(F_nonrandom.*conj(F_nonrandom))))/numel(F_nonrandom);
            %             out_sampled_nonrandom = fftshift(ifftn(ifftshift(F_sampled_nonrandom.*conj(F_sampled_nonrandom))))/numel(F_sampled_nonrandom);
            
            list4(ii,:) = out;
        end
        list4 = list4(:,(size(list4, 2)+1)/2:end);
        list4 = list4./list4(:, 1);
        list5 = list4;
%         list4(list4<0) = NaN;
        %         list44(list44<0) = NaN;
        %                 figure;plot(exp(nanmean(log(list4))));
        %                 hold on;plot(nanmean(list4));plot(acf);hold off;
        
        list4(list4<0) = NaN;
        avg_slopes = diff(nanmean(log(list4)));
        
        temp_list5 = mean(list5);
        temp_list5(temp_list5<0) = NaN;
        avg_exponentials = diff(log(temp_list5));
        clear temp_list5
        
        %         acf = acf(1:sampling:end);
        %         theory_slope = diff(log(acf((length(acf)+1)/2:end)));
        dt = t_s(2) - t_s(1);
        
        list2(i, :) = [slope2D(avg_exponentials(2), dt, k), slope2D(avg_slopes(2), dt, k)];
        
        
        %         figure;plot(t_s(1:end-1), slope2D(avg_exponentials, dt, k));hold on;
        %         plot(t_s(1:end-1), slope2D(avg_slopes, dt, k));
        %         plot(t_s(1:end-1), slope2D(theory_slope, dt, k));
        %
        %                 figure;plot(t(1:end-17), slope2D(avg_exponentials, dt, k));hold on;
        %         plot(t(1:end-17), slope2D(avg_slopes, dt, k));
        %         plot(t(1:end-17), slope2D(theory_slope(1:end-16), dt, k));
        % %
        %         title(['D = ', num2str(D(d)), '; Sampling + Exposure']);
        %         xlabel('time (s)');
        %         ylabel('Diffusion Coefficent (um^2/s)');
        %         xlim([0,0.3]);
        %         ylim([-1.2*D(d),1.2*D(d)]);
        %         legend('Averaging Slopes', 'Averaging Exponentials', 'Theory');
        
        
        
        %
        %         list2(i,:) = [...
        %             (4*k^2*(find_correlation_decay_set_range(t, nanmean(list4)))^-1)^-1,...
        %             (4*k^2*(find_correlation_decay_set_range(t, exp(nanmean(log(list4)))))^-1)^-1,...
        %             (4*k^2*(find_correlation_decay_set_range(t, nanmean(list44)))^-1)^-1,...
        %             (4*k^2*(find_correlation_decay_set_range(t, exp(nanmean(log(list44)))))^-1)^-1];
        
        %         list2(i,:) = [find_correlation_decay(t, out), find_correlation_decay(t, out_smooth)];
        %                     list2(i, :) = [var(signal), var(smooth(signal, exposure_time/dt))];  %Testing exposure time
        % list2(i, :) = [var(signal), var(smoothed_signal), var(signal(1:sampling:end)), var(smoothed_signal(1:sampling:end))];  %Testing exposure time
        % list2(i, :) = [var(signal), var(smoothed_signal), out(1)];
        % list2(i, :) = [var(signal), var(signal(1:sampling(d):end))];
        
        %             list(i,:) = out;
    end
    %     find_correlation_decay(t, mean(list2))
    %     list3(d,:) = [find_correlation_decay(t, mean(list2)), find_correlation_decay(t, mean(list22))];
    list3(d,:) = [nanmean(list2(:,1)), nanmean(list2(:,2)), nanvar(list2(:,1)), nanvar(list2(:,2))];%, mean(list2(:,3)), mean(list2(:,4))];
    %     list3(d,:) = [mean(list2(:,1)), mean(list2(:,2)), mean(list2(:,3)), mean(list2(:,4))];%, mean(list2(:,3)), mean(list2(:,4))];
    
    
end
% hold on
% out = mean(list);

% mean(list(:,(length(out)+1)/2))
% plot(list3(:,1));hold on; plot(list3(:,2));
% figure;plot(D, list3(:,1));hold on; plot(D, list3(:,3));plot(D, list3(:,2));plot(D, list3(:,4));

end

function final_waveform = genarate_random_1D_signal(acf, A)

psd  = abs(fftshift(fftn(acf)));

% seed the random number generator
rng('shuffle');

sample_waveform = ifftn(ifftshift(... % take 3-D inverse FFT of the randomly selected frequency array (ifftshift first to get the center frequency back to the corner)
    randn(size(psd)).*... %(ones(size(K))).*... %  randomly selected frequency points with Gaussian distribution
    sqrt(A.*psd))); % A^3*each frequency point is independent with power (A*B*C*power_spec(a,b,c))& the power spectrum part

% Add real and imaginary parts to get non-symmetric full array:
final_waveform = real(sample_waveform)+imag(sample_waveform);

end

function final_waveform = genarate_nonrandom_1D_signal(acf, A)

psd  = abs(fftshift(fftn(acf)));

sample_waveform = ifftn(ifftshift(... % take 3-D inverse FFT of the randomly selected frequency array (ifftshift first to get the center frequency back to the corner)
    (ones(size(psd))).*... %  randomly selected frequency points with Gaussian distribution
    sqrt(A.*psd))); % A^3*each frequency point is independent with power (A*B*C*power_spec(a,b,c))& the power spectrum part

% Add real and imaginary parts to get non-symmetric full array:
final_waveform = real(sample_waveform)+imag(sample_waveform);

end

function decay_coeff = find_correlation_decay(t, acf)
acf = acf((length(acf)+1)/2:end); %Remove symetric part of acf
acf(acf<0) = NaN; % Remove negative values for calculation log of acf
ln_acf = log(acf); % calc ln(acf) to find slope of this line
diff2_ln_acf = abs(diff(diff(log(acf)))); % 2nd derivative of ln(acf)

% Find linear region of ln(acf) to find slope
% lim1=find(diff2_ln_acf>0.05,1)+1;
lim1 = 20;
lim2=find(isnan(diff2_ln_acf),1);
if isempty(lim2)
    lim2=length(diff2_ln_acf);
end
lim=min(lim1,lim2)-1;

if lim < 2
    lim
    decay_coeff = NaN;
else
    p = polyfit(t(1:lim),ln_acf(1:lim),1);
    
    %         line_fit = polyval(p, t);
    %         plot(t, line_fit)
    %         xlim([0 0.5]);
    %         hold off
    %         pause()
    
    % decay_coeff = -p(1)/(k^2);
    decay_coeff = -p(1);
end
end

function decay_coeff = find_correlation_decay_set_range(t, acf, range)

if mod(length(t),2) == 0
    acf = acf((length(acf))/2 + 1:end); %Remove symetric part of acf
    t = t((length(t))/2 + 1:end);
else
    acf = acf((length(acf)+1)/2:end); %Remove symetric part of acf
    t = t((length(t)+1)/2:end);
end
acf(acf<0) = NaN; % Remove negative values for calculation log of acf
ln_acf = log(acf); % calc ln(acf) to find slope of this line
% diff2_ln_acf = abs(diff(diff(log(acf)))); % 2nd derivative of ln(acf)


% Find linear region of ln(acf) to find slope
% lim1=find(diff2_ln_acf>0.05,1)+1;
lim1 = range(2);
lim2=find(isnan(ln_acf(range(1):end)),1);
if isempty(lim2)
    lim2=lim1;
end
lim=min(lim1,lim2)

if lim-range(1) < 2
    lim
    decay_coeff = NaN;
else
    p = polyfit(t(range(1):lim),ln_acf(range(1):lim),1);
    
    %         line_fit = polyval(p, t);
    %         plot(t(30:100), line_fit(30:100))
    %         xlim([0 0.5]);
    %         hold off
    %         pause()
    
    decay_coeff = -p(1);
    % decay_coeff = NaN;
end
end

function d_list = slope2D(slopes, dt, k)
slopes = slopes/dt;
d_list = (4*k^2*(slopes).^-1).^-1;
end

function filtered_signal = sample_and_exposure_signal(signal, exposure)
% Assuming sampling and exposure window are the same
filtered_signal = zeros(1, floor(length(signal)/exposure));
for i = 1 : floor(length(signal)/exposure)
    filtered_signal(i) = sum(signal(1 + exposure*(i-1):exposure*(i)));
end
end

function filtered_signal = sample_signal(signal, exposure)
% Assuming sampling and exposure window are the same
filtered_signal = zeros(1, floor(length(signal)/exposure));
for i = 1 : floor(length(signal)/exposure)
    filtered_signal(i) = signal(1 + exposure*(i-1));
end
end

function filtered_signal = exposure_signal(signal, exposure)
% Assuming sampling and exposure window are the same
filtered_signal = zeros(1, floor(length(signal)/exposure));
for i = 1 : length(signal)-exposure
    filtered_signal(i) = sum(signal(i:i+exposure));
end
end
