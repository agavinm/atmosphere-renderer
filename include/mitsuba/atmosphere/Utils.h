#ifndef MITSUBA_UTILS_H
#define MITSUBA_UTILS_H

#include <array>

namespace Utils {

    /*template <typename Float, typename UInt32, typename Mask, std::size_t SIZE>
    UInt32 get_low(const std::array<float, SIZE> &keys, const Float &x) {
        UInt32 index = 0;
        UInt32 begining = 0;
        UInt32 ending = keys.size() - 1;
        Mask not_found = x >= Float(keys[0]) && x <= Float(keys.back());

        while (enoki::any(not_found)) {
            index = enoki::select(
                    not_found,
                    begining + (ending - begining) / UInt32(2),
                    index
            );

            const Mask msk_if = enoki::gather<Float>(keys.data(), index) <= x && enoki::gather<Float>(keys.data(), index + 1) >= x;
            not_found = enoki::select(
                    not_found && msk_if,
                    !not_found, // Found
                    not_found // Not found yet
            );

            // Dicotomic search
            const Mask msk_if2 = x > enoki::gather<Float>(keys.data(), index);
            begining = enoki::select(
                    !msk_if && msk_if2,
                    index,
                    begining
            );
            ending = enoki::select(
                    !msk_if && !msk_if2,
                    index,
                    ending
            );
        }

        return index;
    }

    template<std::size_t SIZE>
    size_t get_low(const std::array<float, SIZE> &keys, const float &x) {
        size_t index = 0;
        size_t begining = 0;
        size_t ending = keys.size() - 1;

        while (true) {
            index = begining + (ending - begining) / 2;
            if (keys[index] <= x && keys[index+1] >= x) {
                return index; // That's the index that i'm looking for!
            }

            else {
                // Dicotomic search
                if (x > keys[index]) {
                    begining = index;
                }
                else {
                    ending = index;
                }
            }
        }
    }

    template <typename Float, typename UInt32, typename Mask, std::size_t SIZE>
    Float interpolate(const std::array<float, SIZE> &keys, const std::array<float, SIZE> &values, const Float &x) {
        const UInt32 idx = get_low<Float, UInt32, Mask>(keys, x);

        const Float keys_idx = enoki::gather<Float>(keys.data(), idx);
        const Float weight = (x - keys_idx) / (enoki::gather<Float>(keys.data(), idx + 1) - keys_idx);
        const Float finalResult = enoki::select(
                x > Float(keys.back()),
                Float(values.back()),
                enoki::select(
                        x < Float(keys[0]),
                        Float(values[0]),
                        (Float(1.) - weight) * enoki::gather<Float>(values.data(), idx) + weight * enoki::gather<Float>(values.data(), idx + 1) // Is valid
                )
        );

        return finalResult;
    }

    template <typename Float = float, typename, typename, std::size_t SIZE>
    [[nodiscard]] float interpolate(const std::array<float, SIZE> &keys, const std::array<float, SIZE> &values, const float &x) {
        if (x < keys[0])
            return values[0];

        if (x >= keys.back())
            return values.back();

        const size_t idx = get_low(keys, x);
        const float weight = (x - keys[idx]) / (keys[idx + 1] - keys[idx]);
        const float finalResult = (1.f - weight) * values[idx] + weight * values[idx + 1];

        return finalResult;
    }*/

    template <typename Float>
    Float fast_interpolate(const Float &wl_lower_limit, const Float &wl_upper_limit, const Float &wl_x, const Float &lower_limit, const Float &upper_limit) {
        const Float weight = (wl_x - wl_lower_limit) / (wl_upper_limit - wl_lower_limit);
        return ((Float(1.) - weight) * lower_limit ) +  (weight * upper_limit);
    }

    // Return max value
    /*template<std::size_t SIZE>
    float max(const std::array<float, SIZE> &values) {
        if (values.empty())
            return 0;

        float finalResult = values[0];
        for (size_t i = 1; i < values.size(); i++) {
            if (values[i] > finalResult)
                finalResult = values[i];
        }

        return finalResult;
    }*/

    template <typename Float, typename ScalarFloat>
    ScalarFloat get_first(const Float &x) {
        return enoki::slice(x, 0);
    }

    template <typename Float = float, typename>
    [[nodiscard]] float get_first(const float &x) {
        return x;
    }

    template <typename Float = double, typename>
    [[nodiscard]] double get_first(const double &x) {
        return x;
    }

    template <typename Float, typename UInt32, typename Mask, std::size_t SIZE>
    Float get_(const Float &key, const std::array<float, SIZE> &values) {
        const Mask msk1 = key < Float(0), msk2 = key >= Float(values.size());
        const UInt32 index = enoki::select(
                msk1,
                UInt32(0),
                enoki::select(
                        msk2,
                        UInt32(values.size() - 1),
                        UInt32(key)
                        )
        );

        const UInt32 upper_index = enoki::select(
                index < UInt32(values.size() - 1),
                index + UInt32(1),
                UInt32(0)
        );

        return enoki::select(
                msk1 || msk2 || key == Float(index),
                enoki::gather<Float>(values.data(), index),
                Utils::fast_interpolate(Float(index),
                                        Float(upper_index),
                                        key,
                                        enoki::gather<Float>(values.data(), index),
                                        enoki::gather<Float>(values.data(), upper_index))
        );
    }

    template <typename Float = float, typename, typename, std::size_t SIZE>
    [[nodiscard]] float get_(const float &key, const std::array<float, SIZE> &values) {
        if (key < 0.f)
            return values[0];
        else if (key >= float(values.size()))
            return values.back();

        const int index = int(key);
        const auto index_float = float(index);
        if (key == index_float)
            return values[index];

        return Utils::fast_interpolate<float>(index_float, index_float + 1.f, key, values[index], values[index + 1]);
    }
}; //Utils

#endif // MITSUBA_UTILS_H
