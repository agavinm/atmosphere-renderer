#ifndef MITSUBA_UTILS_H
#define MITSUBA_UTILS_H

#include <vector>

namespace Utils {

template <typename Float, typename UInt32, typename Mask>
UInt32 get_low(const std::vector<float> &keys, const Float &x) {
    UInt32 index = 0;
    UInt32 begining = 0;
    UInt32 ending = keys.size() - 1;
    Mask not_found = x >= Float(keys[0]) && x <= Float(keys[keys.size() - 1]);

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

size_t get_low(const std::vector<float> &keys, const float &x) {
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

template <typename Float, typename UInt32, typename Mask>
Float interpolate(const std::vector<float> &keys, const std::vector<float> &values, const Float &x) {
    const UInt32 idx = get_low<Float, UInt32, Mask>(keys, x);

    const Float keys_idx = enoki::gather<Float>(keys.data(), idx);
    const Float weight = (x - keys_idx) / (enoki::gather<Float>(keys.data(), idx + 1) - keys_idx);
    const Float finalResult = enoki::select(
            x > Float(keys[keys.size() - 1]),
            Float(values[values.size() - 1]),
            enoki::select(
                    x < Float(keys[0]),
                    Float(values[0]),
                    (Float(1.) - weight) * enoki::gather<Float>(values.data(), idx) + weight * enoki::gather<Float>(values.data(), idx + 1) // Is valid
                    )
            );

    return finalResult;
}

template <typename Float = float, typename, typename>
float interpolate(const std::vector<float> &keys, const std::vector<float> &values, const float &x) {
    if (x < keys[0])
        return values[0];

    if (x >= keys.at(keys.size()-1))
        return values[keys.size() - 1];

    const size_t idx = get_low(keys, x);
    const float weight = (x - keys[idx]) / (keys[idx + 1] - keys[idx]);
    const float finalResult = (1.f - weight) * values[idx] + weight * values[idx + 1];

    return finalResult;
}

template <typename Float>
Float fast_interpolate(const Float &wl_lower_limit, const Float &wl_upper_limit, const Float &wl_x, const Float &lower_limit, const Float &upper_limit) {
    const Float weight = (wl_x - wl_lower_limit) / (wl_upper_limit - wl_lower_limit);
    return ((Float(1.) - weight)*lower_limit ) +  (weight*upper_limit);
}

// Return max value
float max(const std::vector<float> &values) {
    if (values.empty())
        return 0;

    float finalResult = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        if (values[i] > finalResult)
            finalResult = values[i];
    }

    return finalResult;
}

}; //Utils

#endif // MITSUBA_UTILS_H
