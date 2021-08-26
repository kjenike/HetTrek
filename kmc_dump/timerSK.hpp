#ifndef TIMER_SK_HPP
#define TIMER_SK_HPP
#include <string>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <cassert>

class TimerSK {
    private:
        std::chrono::high_resolution_clock::time_point start;
    public:
        inline TimerSK() {
            reset();
        }
        inline void reset() {
            start = std::chrono::high_resolution_clock::now();
        }
        inline double get() {
            return (std::chrono::duration_cast< std::chrono::duration<double> > (std::chrono::high_resolution_clock::now() - start).count()) * 1000.0;
        }
        inline double lap() {
            double ret = get();
            reset();
            return ret;
        }
};
#endif
