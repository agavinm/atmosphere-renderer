/* 
* Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom
* the Software is furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
* OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef MITSUBA_UTILS_H
#define MITSUBA_UTILS_H

#include <vector>

namespace Utils {

template <typename Float>
unsigned int get_low(const std::vector<Float> keys, const Float &x) {
    unsigned int index = 0;
    Float begining = 0;
    Float ending = keys.size() - 1;

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

template <typename Float>
Float interpolate(const std::vector<Float> keys, const Float* values, const Float &x) {
    if (x < keys[0])
        return values[0];

    Float last;
    last = keys.at(keys.size()-1);

    if (x >= last) {
        //Float value = *keys.end();
        //printf("Keys last value : %f\n", value);
        return values[keys.size() - 1];
    }

    unsigned int idx = get_low(keys, x);
    Float weight = (x - keys[idx]) / (keys[idx + 1] - keys[idx]);
    Float finalResult = (1. - weight)*values[idx] + weight*values[idx + 1];
    
    return finalResult;
}

template <typename Float>
Float fast_interpolate(Float wl_lower_limit, Float wl_upper_limit, Float wl_x, Float lower_limit, Float upper_limit) {
    Float weight = (wl_x - wl_lower_limit) / (wl_upper_limit - wl_lower_limit);
    Float finalResult = ((1. - weight)*lower_limit ) +  (weight*upper_limit);

    return finalResult;
}

}; //Utils

#endif // MITSUBA_UTILS_H
