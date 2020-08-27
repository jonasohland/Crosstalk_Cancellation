//
//  CrossTalkCancellation.hpp
//  CAPlayThrough
//
//  Created by MINGFENWANG on 2018/9/14.
//

#ifndef CrossTalkCancellation_hpp
#define CrossTalkCancellation_hpp

#include <vector>

template <typename SampleType>
class CrossTalkCancellation {
public:
    CrossTalkCancellation() {
        sr = 44100;
        spkr_to_spkr = 0.3048;
        lstnr_to_spkr = 0.5588;
        ear_to_ear = 0.215;
        headshadow_filter_coefficients(compute_geometry(), ear_to_ear/2.0);

        // init left queue
        for (int i = 0; i < 4; ++i) {
            left_queue.push_back(std::make_unique<std::vector<SampleType>>());
            left_queue[i]->resize(512, 0.0);
            right_queue.push_back(std::make_unique<std::vector<SampleType>>());
            right_queue[i]->resize(512, 0.0);
        }
    }
    
    void enable()
    {
        enabled = true;
    }
    
    void disable()
    {
        enabled = false;
    }
    
    void set_blocksize(int size) {
        left_queue.clear();
        right_queue.clear();
        for (int i = 0; i < 4; ++i) {
            left_queue.push_back(std::make_unique<std::vector<SampleType>>());
            left_queue[i]->resize(size, 0.0);
            right_queue.push_back(std::make_unique<std::vector<SampleType>>());
            right_queue[i]->resize(size, 0.0);
        }
    }
    
    void recalc_filters()
    {
        headshadow_filter_coefficients(compute_geometry(), ear_to_ear/2.0);
    }
    
    void set_samplerate(const int _sr) {
        sr = _sr;
    }
    
    void set_speaker_distance(const SampleType dist) {
        spkr_to_spkr = dist;
        
    }
    
    void set_listener_distance(const SampleType dist) {
        lstnr_to_spkr = dist;
    }
    
    void set_ear_distance(const SampleType dist) {
        ear_to_ear = dist;
    }
    
    void process_stereo_channel(SampleType* left, SampleType* right, std::size_t samples) {
        
        if (!enabled)
            return;
        
        // move old data pieces forward
        left_queue[0] = std::move(left_queue[1]);
        left_queue[1] = std::move(left_queue[2]);
        left_queue[2] = std::move(left_queue[3]);
        // queue new data piece
        left_queue[3] = std::make_unique<std::vector<SampleType>>(left, left + samples);

        // move old data pieces forward
        right_queue[0] = std::move(right_queue[1]);
        right_queue[1] = std::move(right_queue[2]);
        right_queue[2] = std::move(right_queue[3]);
        // queue new data piece
        right_queue[3] = std::make_unique<std::vector<SampleType>>(right, right + samples);

        // concatenate data pieces to one big data piece
        std::vector<SampleType> work_left;
        work_left.insert(work_left.end(), left_queue[0]->begin(), left_queue[0]->end());
        work_left.insert(work_left.end(), left_queue[1]->begin(), left_queue[1]->end());
        work_left.insert(work_left.end(), left_queue[2]->begin(), left_queue[2]->end());
        work_left.insert(work_left.end(), left_queue[3]->begin(), left_queue[3]->end());

        std::vector<SampleType> work_right;
        work_right.insert(work_right.end(), right_queue[0]->begin(), right_queue[0]->end());
        work_right.insert(work_right.end(), right_queue[1]->begin(), right_queue[1]->end());
        work_right.insert(work_right.end(), right_queue[2]->begin(), right_queue[2]->end());
        work_right.insert(work_right.end(), right_queue[3]->begin(), right_queue[3]->end());

        // calculate crosstalk cancellation for left channel
        std::vector<SampleType> l_left, l_right;
        l_left.resize(work_left.size(), 0.0);
        l_right.resize(work_left.size(), 0.0);
        cancel_crosstalk(work_left, l_left, l_right);
        
        // calculate crosstalk cancellation for right channel
        std::vector<SampleType> r_right, r_left;
        r_right.resize(work_right.size(), 0.0);
        r_left.resize(work_right.size(), 0.0);
        cancel_crosstalk(work_right, r_right, r_left);
        
        // accumulate crosstalk cancellation to left channel
        add_to_signal(work_left, l_left);
        add_to_signal(work_left, r_left);
        
        // accumulate crosstalk cancellation to right channel
        add_to_signal(work_right, r_right);
        add_to_signal(work_right, l_right);
        
        // store the third data piece as the result
        
        std::size_t data_size = (left_queue[2]->size() < samples) ? left_queue[2]->size() : samples;
       
        memcpy(left,
               work_left.data() + left_queue[0]->size() + left_queue[1]->size(),
               data_size * sizeof(SampleType));
        
        memcpy(right,
               work_right.data() + right_queue[0]->size() + right_queue[1]->size(),
               data_size * sizeof(SampleType));
    }
private:
    /******************************************************************************
     IIR(Infinite impulse response) filter
     Where:
     b[i] are the feedforward filter coefficients.
     a[i] are the feedback filter coefficients.
     *****************************************************************************/
    struct HeadShadow {
        SampleType b[2];
        SampleType a[2];
    };
    
    
    void cancel_crosstalk(const std::vector<SampleType>& signal,
                          std::vector<SampleType>& ipsilateral,
                          std::vector<SampleType>& contralateral) {
        SampleType c = 343.2;
        SampleType delta_d = abs(d2 - d1);
        SampleType time_delay = delta_d / c;
        SampleType attenuation = d1 / d2;
        // Reference max amplitude
        SampleType ref = max_of_abs(signal);
        recursive_cancel(signal, ipsilateral, contralateral, ref, time_delay, attenuation);
    }
    
    inline void add_to_signal(std::vector<SampleType>& signal,
                       const std::vector<SampleType>& signal2) {
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            signal[i] += signal2[i];
        }
    }
    
    /******************************************************************************
                ----------------        ----------------
               |   M samples    |      |delta fractional|
     y[n] ---> | integer delay  | ---> |      delay     | ---> y[n - M - delta]
               |                |      |                |
                ----------------        ----------------
     
     The M samples integer delay is implemented as a simple delay line:
     y[n] = x[n - M]
     
     The delta fractional delay is implemented as:
     y[n] = u * x[n] + x[n - 1] - u * y[n - 1] where u = (1 - delta) / (1 + delta)
     
     The formula has initial condition issues, so I use naive method instead.
     *****************************************************************************/
    inline void fractional_delay(std::vector<SampleType>& signal,
                          SampleType time) {
        std::vector<SampleType> temp_signal(signal.size());
        SampleType m = time * sr;
        int m_int = floor(m);
        SampleType m_frac = m - m_int;
        for (int i=0; i<static_cast<int>(temp_signal.size()); i++) {
            int index_low = i - (m_int + 1);
            int index_high = i - m_int;
            SampleType low = (index_low >= 0 ? signal[index_low] : 0.0);
            SampleType high = (index_high >= 0 ? signal[index_high] : 0.0);
            temp_signal[i] = high + (low - high) * m_frac;
        }
        signal = temp_signal;
    }
    
    /******************************************************************************
     A linear filter that achieves zero phase delay by applying an IIR filter to a
     signal twice, once forwards and once backwards. The order of the filter is
     twice the original filter order.
     y[t] = b[0] * x[t] + b[1] * x[t-1] + ... - a[1] * y[t-1] - a[2] * y[t-2] - ...
     To simplify the caculation, I let y[0] = x[0], though this method has initial
     condition issue, I use padding data at both head and tail to avoid the issue.
     *****************************************************************************/
    inline void filtfilt(std::vector<SampleType>& signal) {
        if (signal.size() > 1) {
            std::vector<SampleType> temp_signal;
            // from left to right
            temp_signal = signal;
            for (int i=1; i<static_cast<int>(signal.size()); i++) {
                temp_signal[i] = headshadow.b[0] * signal[i] +
                headshadow.b[1] * signal[i-1] -
                headshadow.a[1] * temp_signal[i-1];
            }
            
            // from right to left
            signal = temp_signal;
            for (int i=static_cast<int>(signal.size())-2; i>=0; i--) {
                signal[i] = headshadow.b[0] * temp_signal[i] +
                headshadow.b[1] * temp_signal[i+1] -
                headshadow.a[1] * signal[i+1];
            }
        }
    }
    
    inline void attenuate(std::vector<SampleType>& signal, const SampleType attenuation) {
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            signal[i] *= attenuation;
        }
    }
    
    /******************************************************************************
     To speed up calculation and optimize memory, I use loop instead of recursive,
     accumulate crosstalk cancellation to each channel in turn.
     *****************************************************************************/
    void recursive_cancel(const std::vector<SampleType>& signal,
                          std::vector<SampleType>& ipsilateral,
                          std::vector<SampleType>& contralateral,
                          const SampleType ref,
                          const SampleType time,
                          const SampleType attenuation,
                          const SampleType threshold_db = -70.0) {
        std::vector<SampleType> temp_signal = signal;
        bool ping_pong = false;
        SampleType db;
        do {
            // time delay by linear interpolation
            fractional_delay(temp_signal, time);

            // invert the delayed signal
            invert(temp_signal);

            // apply headshadow filter (lowpass based on theta)
            filtfilt(temp_signal);

            // attenuate the low pass filtered delayed signal
            attenuate(temp_signal, attenuation);

            // accumulate signal to either ipsilateral or contralateral
            add_to_signal(ping_pong ? ipsilateral : contralateral,
                          temp_signal);

            // Recurse until rms db is below threshold
            db = 20 * log10(max_of_abs(temp_signal) / ref);

            // flip left and right
            ping_pong = !ping_pong;
        } while (db >= threshold_db);
    }
    
    inline SampleType max_of_abs(const std::vector<SampleType>& signal) {
        SampleType max_abs_x = 0.0;
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            SampleType abs_x = std::abs(signal[i]);
            if (abs_x > max_abs_x) {
                max_abs_x = abs_x;
            }
        }
        return max_abs_x;
    }
    
    void headshadow_filter_coefficients(const SampleType _theta,
                                        const SampleType r) {
        SampleType theta = _theta + M_PI / 2.0;
        SampleType theta0 = 2.618;
        SampleType alpha_min = 0.5;
        SampleType c = 343.2;
        SampleType w0 = c / r;
        SampleType alpha = 1 + alpha_min / 2.0 + (1.0 - alpha_min / 2.0) * cos(theta * M_PI / theta0);
        
        headshadow.b[0] = (alpha + w0 / sr) / (1 + w0 / sr);
        headshadow.b[1] = (-alpha + w0 / sr) / (1 + w0 / sr);
        headshadow.a[0] = 1;
        headshadow.a[1] = -(1 - w0 / sr) / (1 + w0 / sr);
    }
    
    SampleType compute_geometry() {
        SampleType S = spkr_to_spkr / 2.0;
        SampleType L = lstnr_to_spkr;
        SampleType r = ear_to_ear / 2.0;
        SampleType theta = acos(S / (sqrt(L * L + S * S)));
        SampleType delta_d = r * (M_PI - 2.0 * theta);
        d1 = sqrt(L * L + (S - r) * (S - r));
        d2 = d1 + delta_d;
        
        // angle from center of head to speaker (used for computing headshadow)
        theta = atan(S / L);
        
        return theta;
    }
    
    inline void invert(std::vector<SampleType>& signal) {
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            signal[i] = -signal[i];
        }
    }
    
    bool enabled = true;
    
    SampleType d1;
    SampleType d2;
    HeadShadow headshadow;
    SampleType spkr_to_spkr;
    SampleType lstnr_to_spkr;
    SampleType ear_to_ear;
    int sr;
    std::vector< std::unique_ptr<std::vector<SampleType>> > left_queue;
    std::vector< std::unique_ptr<std::vector<SampleType>> > right_queue;
};

#endif /* CrossTalkCancellation_hpp */
