#ifndef _AEROSOL_POLAR_ANTARTIC
#define _AEROSOL_POLAR_ANTARTIC

#include <array>
#include <mitsuba/atmosphere/GlobalAtmosphericAerosol.h>

namespace polarAntarticAerosol {

	// row[0] = wavelength
	// row[1] = cross section absorption coefficient
	// row[2] = cross section scattering coefficient
    const static std::array<std::array<float, 1001>, 3> tabulatedValues =
	{ std::array<float, 1001>{ 100,100.23,100.46,100.69,100.93,101.16,101.39,101.62,101.86,102.09,102.33,102.57,102.8,103.04,103.28,103.51,103.75,103.99,104.23,104.47,104.71,104.95,105.2,105.44,105.68,105.93,106.17,106.41,106.66,106.91,107.15,107.4,107.65,107.89,108.14,108.39,108.64,108.89,109.14,109.4,109.65,109.9,110.15,110.41,110.66,110.92,111.17,111.43,111.69,111.94,112.2,112.46,112.72,112.98,113.24,113.5,113.76,114.03,114.29,114.55,114.82,115.08,115.35,115.61,115.88,116.14,116.41,116.68,116.95,117.22,117.49,117.76,118.03,118.3,118.58,118.85,119.12,119.4,119.67,119.95,120.23,120.5,120.78,121.06,121.34,121.62,121.9,122.18,122.46,122.74,123.03,123.31,123.59,123.88,124.17,124.45,124.74,125.03,125.31,125.6,125.89,126.18,126.47,126.77,127.06,127.35,127.64,127.94,128.23,128.53,128.82,129.12,129.42,129.72,130.02,130.32,130.62,130.92,131.22,131.52,131.83,132.13,132.43,132.74,133.05,133.35,133.66,133.97,134.28,134.59,134.9,135.21,135.52,135.83,136.14,136.46,136.77,137.09,137.4,137.72,138.04,138.36,138.68,139,139.32,139.64,139.96,140.28,140.6,140.93,141.25,141.58,141.91,142.23,142.56,142.89,143.22,143.55,143.88,144.21,144.54,144.88,145.21,145.55,145.88,146.22,146.55,146.89,147.23,147.57,147.91,148.25,148.59,148.94,149.28,149.62,149.97,150.31,150.66,151.01,151.36,151.71,152.05,152.41,152.76,153.11,153.46,153.82,154.17,154.53,154.88,155.24,155.6,155.96,156.31,156.68,157.04,157.4,157.76,158.12,158.49,158.85,159.22,159.59,159.96,160.32,160.69,161.06,161.44,161.81,162.18,162.55,162.93,163.31,163.68,164.06,164.44,164.82,165.2,165.58,165.96,166.34,166.72,167.11,167.49,167.88,168.27,168.66,169.04,169.43,169.82,170.22,170.61,171,171.4,171.79,172.19,172.58,172.98,173.38,173.78,174.18,174.58,174.98,175.39,175.79,176.2,176.6,177.01,177.42,177.83,178.24,178.65,179.06,179.47,179.89,180.3,180.72,181.13,181.55,181.97,182.39,182.81,183.23,183.65,184.08,184.5,184.93,185.35,185.78,186.21,186.64,187.07,187.5,187.93,188.36,188.8,189.23,189.67,190.11,190.55,190.99,191.43,191.87,192.31,192.75,193.2,193.64,194.09,194.54,194.98,195.43,195.88,196.34,196.79,197.24,197.7,198.15,198.61,199.07,199.53,199.99,200.45,200.91,201.37,201.84,202.3,202.77,203.24,203.7,204.17,204.64,205.12,205.59,206.06,206.54,207.01,207.49,207.97,208.45,208.93,209.41,209.89,210.38,210.86,211.35,211.84,212.32,212.81,213.3,213.8,214.29,214.78,215.28,215.77,216.27,216.77,217.27,217.77,218.27,218.78,219.28,219.79,220.29,220.8,221.31,221.82,222.33,222.84,223.36,223.87,224.39,224.91,225.42,225.94,226.46,226.99,227.51,228.03,228.56,229.09,229.61,230.14,230.67,231.21,231.74,232.27,232.81,233.35,233.88,234.42,234.96,235.5,236.05,236.59,237.14,237.68,238.23,238.78,239.33,239.88,240.44,240.99,241.55,242.1,242.66,243.22,243.78,244.34,244.91,245.47,246.04,246.6,247.17,247.74,248.31,248.89,249.46,250.03,250.61,251.19,251.77,252.35,252.93,253.51,254.1,254.68,255.27,255.86,256.45,257.04,257.63,258.23,258.82,259.42,260.02,260.62,261.22,261.82,262.42,263.03,263.63,264.24,264.85,265.46,266.07,266.69,267.3,267.92,268.53,269.15,269.77,270.4,271.02,271.64,272.27,272.9,273.53,274.16,274.79,275.42,276.06,276.69,277.33,277.97,278.61,279.25,279.9,280.54,281.19,281.84,282.49,283.14,283.79,284.45,285.1,285.76,286.42,287.08,287.74,288.4,289.07,289.73,290.4,291.07,291.74,292.42,293.09,293.76,294.44,295.12,295.8,296.48,297.17,297.85,298.54,299.23,299.92,300.61,301.3,302,302.69,303.39,304.09,304.79,305.49,306.2,306.9,307.61,308.32,309.03,309.74,310.46,311.17,311.89,312.61,313.33,314.05,314.77,315.5,316.23,316.96,317.69,318.42,319.15,319.89,320.63,321.37,322.11,322.85,323.59,324.34,325.09,325.84,326.59,327.34,328.1,328.85,329.61,330.37,331.13,331.89,332.66,333.43,334.19,334.97,335.74,336.51,337.29,338.06,338.84,339.63,340.41,341.19,341.98,342.77,343.56,344.35,345.14,345.94,346.74,347.54,348.34,349.14,349.95,350.75,351.56,352.37,353.18,354,354.81,355.63,356.45,357.27,358.1,358.92,359.75,360.58,361.41,362.24,363.08,363.92,364.75,365.59,366.44,367.28,368.13,368.98,369.83,370.68,371.54,372.39,373.25,374.11,374.97,375.84,376.7,377.57,378.44,379.32,380.19,381.07,381.94,382.82,383.71,384.59,385.48,386.37,387.26,388.15,389.05,389.94,390.84,391.74,392.64,393.55,394.46,395.37,396.28,397.19,398.11,399.02,399.94,400.87,401.79,402.72,403.65,404.58,405.51,406.44,407.38,408.32,409.26,410.2,411.15,412.1,413.05,414,414.95,415.91,416.87,417.83,418.79,419.76,420.73,421.7,422.67,423.64,424.62,425.6,426.58,427.56,428.55,429.54,430.53,431.52,432.51,433.51,434.51,435.51,436.52,437.52,438.53,439.54,440.55,441.57,442.59,443.61,444.63,445.66,446.68,447.71,448.75,449.78,450.82,451.86,452.9,453.94,454.99,456.04,457.09,458.14,459.2,460.26,461.32,462.38,463.45,464.52,465.59,466.66,467.74,468.81,469.89,470.98,472.06,473.15,474.24,475.34,476.43,477.53,478.63,479.73,480.84,481.95,483.06,484.17,485.29,486.41,487.53,488.65,489.78,490.91,492.04,493.17,494.31,495.45,496.59,497.74,498.88,500.03,501.19,502.34,503.5,504.66,505.82,506.99,508.16,509.33,510.51,511.68,512.86,514.04,515.23,516.42,517.61,518.8,520,521.19,522.4,523.6,524.81,526.02,527.23,528.45,529.66,530.88,532.11,533.33,534.56,535.8,537.03,538.27,539.51,540.75,542,543.25,544.5,545.76,547.02,548.28,549.54,550.81,552.08,553.35,554.63,555.9,557.19,558.47,559.76,561.05,562.34,563.64,564.94,566.24,567.54,568.85,570.16,571.48,572.8,574.12,575.44,576.77,578.1,579.43,580.76,582.1,583.45,584.79,586.14,587.49,588.84,590.2,591.56,592.93,594.29,595.66,597.04,598.41,599.79,601.17,602.56,603.95,605.34,606.74,608.13,609.54,610.94,612.35,613.76,615.18,616.6,618.02,619.44,620.87,622.3,623.73,625.17,626.61,628.06,629.51,630.96,632.41,633.87,635.33,636.8,638.26,639.73,641.21,642.69,644.17,645.65,647.14,648.63,650.13,651.63,653.13,654.64,656.15,657.66,659.17,660.69,662.22,663.74,665.27,666.81,668.34,669.88,671.43,672.98,674.53,676.08,677.64,679.2,680.77,682.34,683.91,685.49,687.07,688.65,690.24,691.83,693.43,695.02,696.63,698.23,699.84,701.46,703.07,704.69,706.32,707.95,709.58,711.21,712.85,714.5,716.14,717.79,719.45,721.11,722.77,724.44,726.11,727.78,729.46,731.14,732.82,734.51,736.21,737.9,739.61,741.31,743.02,744.73,746.45,748.17,749.89,751.62,753.36,755.09,756.83,758.58,760.33,762.08,763.84,765.6,767.36,769.13,770.9,772.68,774.46,776.25,778.04,779.83,781.63,783.43,785.24,787.05,788.86,790.68,792.5,794.33,796.16,797.99,799.83,801.68,803.53,805.38,807.24,809.1,810.96,812.83,814.7,816.58,818.46,820.35,822.24,824.14,826.04,827.94,829.85,831.76,833.68,835.6,837.53,839.46,841.4,843.33,845.28,847.23,849.18,851.14,853.1,855.07,857.04,859.01,860.99,862.98,864.97,866.96,868.96,870.96,872.97,874.98,877,879.02,881.05,883.08,885.12,887.16,889.2,891.25,893.31,895.36,897.43,899.5,901.57,903.65,905.73,907.82,909.91,912.01,914.11,916.22,918.33,920.45,922.57,924.7,926.83,928.97,931.11,933.25,935.41,937.56,939.72,941.89,944.06,946.24,948.42,950.6,952.8,954.99,957.19,959.4,961.61,963.83,966.05,968.28,970.51,972.75,974.99,977.24,979.49,981.75,984.01,986.28,988.55,990.83,993.12,995.41,997.7,1000 },
      std::array<float, 1001>{ 1.1964e-16,1.2e-16,1.2014e-16,1.199e-16,1.194e-16,1.1888e-16,1.1848e-16,1.1841e-16,1.1872e-16,1.1924e-16,1.198e-16,1.2013e-16,1.2019e-16,1.1992e-16,1.194e-16,1.1887e-16,1.1851e-16,1.1845e-16,1.187e-16,1.1925e-16,1.1978e-16,1.2016e-16,1.2028e-16,1.2003e-16,1.1959e-16,1.1899e-16,1.1857e-16,1.1839e-16,1.1856e-16,1.1903e-16,1.1969e-16,1.2011e-16,1.2034e-16,1.2019e-16,1.1976e-16,1.1917e-16,1.1865e-16,1.1834e-16,1.1844e-16,1.1878e-16,1.1934e-16,1.1989e-16,1.2028e-16,1.2036e-16,1.201e-16,1.196e-16,1.1901e-16,1.1851e-16,1.1827e-16,1.1841e-16,1.1882e-16,1.1939e-16,1.1996e-16,1.2032e-16,1.2039e-16,1.2008e-16,1.1955e-16,1.1897e-16,1.1853e-16,1.183e-16,1.1839e-16,1.1885e-16,1.1936e-16,1.1996e-16,1.2034e-16,1.2039e-16,1.2017e-16,1.1969e-16,1.1907e-16,1.1856e-16,1.1829e-16,1.1831e-16,1.1861e-16,1.1924e-16,1.1981e-16,1.2032e-16,1.2058e-16,1.2045e-16,1.2005e-16,1.1947e-16,1.1889e-16,1.1838e-16,1.1826e-16,1.1842e-16,1.1886e-16,1.1944e-16,1.2003e-16,1.2049e-16,1.2065e-16,1.205e-16,1.201e-16,1.1949e-16,1.1889e-16,1.1848e-16,1.1833e-16,1.1849e-16,1.1893e-16,1.1956e-16,1.2017e-16,1.2059e-16,1.2087e-16,1.207e-16,1.2023e-16,1.1969e-16,1.191e-16,1.1862e-16,1.1837e-16,1.1848e-16,1.1893e-16,1.1947e-16,1.2014e-16,1.207e-16,1.21e-16,1.2102e-16,1.2065e-16,1.2016e-16,1.195e-16,1.1898e-16,1.1857e-16,1.1849e-16,1.1871e-16,1.1917e-16,1.1978e-16,1.2045e-16,1.2092e-16,1.2118e-16,1.2115e-16,1.2081e-16,1.2024e-16,1.1962e-16,1.1903e-16,1.1862e-16,1.1852e-16,1.1874e-16,1.1925e-16,1.1984e-16,1.2048e-16,1.2102e-16,1.2131e-16,1.213e-16,1.2103e-16,1.2049e-16,1.1993e-16,1.193e-16,1.1881e-16,1.1866e-16,1.1869e-16,1.1904e-16,1.1965e-16,1.2029e-16,1.209e-16,1.2133e-16,1.2147e-16,1.2137e-16,1.2097e-16,1.2036e-16,1.1973e-16,1.1911e-16,1.1869e-16,1.1853e-16,1.187e-16,1.1912e-16,1.1972e-16,1.204e-16,1.2103e-16,1.2145e-16,1.2164e-16,1.2152e-16,1.2111e-16,1.2054e-16,1.1985e-16,1.1924e-16,1.1879e-16,1.1851e-16,1.1853e-16,1.1886e-16,1.1943e-16,1.2005e-16,1.2078e-16,1.2133e-16,1.2161e-16,1.2161e-16,1.2134e-16,1.2086e-16,1.2019e-16,1.1944e-16,1.1885e-16,1.1852e-16,1.1848e-16,1.1877e-16,1.1923e-16,1.1989e-16,1.2062e-16,1.2121e-16,1.2157e-16,1.2169e-16,1.2152e-16,1.2109e-16,1.2048e-16,1.1978e-16,1.1911e-16,1.1861e-16,1.1839e-16,1.1843e-16,1.1878e-16,1.1934e-16,1.2001e-16,1.2071e-16,1.2131e-16,1.2169e-16,1.2182e-16,1.2166e-16,1.2126e-16,1.2071e-16,1.1998e-16,1.1928e-16,1.1873e-16,1.1836e-16,1.1829e-16,1.1849e-16,1.1899e-16,1.1971e-16,1.204e-16,1.211e-16,1.2161e-16,1.2191e-16,1.2198e-16,1.2167e-16,1.2119e-16,1.2048e-16,1.1978e-16,1.1908e-16,1.1864e-16,1.1835e-16,1.1842e-16,1.1873e-16,1.1921e-16,1.1989e-16,1.2056e-16,1.2126e-16,1.2175e-16,1.2203e-16,1.2208e-16,1.2185e-16,1.2139e-16,1.2075e-16,1.2012e-16,1.1937e-16,1.1886e-16,1.1847e-16,1.1841e-16,1.1855e-16,1.1895e-16,1.1956e-16,1.2032e-16,1.2103e-16,1.2175e-16,1.2226e-16,1.2248e-16,1.2245e-16,1.2213e-16,1.216e-16,1.2091e-16,1.2018e-16,1.1947e-16,1.1888e-16,1.1852e-16,1.1839e-16,1.1851e-16,1.1898e-16,1.1957e-16,1.2034e-16,1.2114e-16,1.2181e-16,1.2236e-16,1.2267e-16,1.227e-16,1.2252e-16,1.2215e-16,1.2156e-16,1.2082e-16,1.2011e-16,1.1948e-16,1.1897e-16,1.1866e-16,1.1869e-16,1.1888e-16,1.1932e-16,1.1987e-16,1.2062e-16,1.2137e-16,1.2205e-16,1.2258e-16,1.229e-16,1.2305e-16,1.229e-16,1.225e-16,1.2205e-16,1.2131e-16,1.2059e-16,1.1983e-16,1.1928e-16,1.1882e-16,1.1861e-16,1.1877e-16,1.1915e-16,1.197e-16,1.2045e-16,1.213e-16,1.2206e-16,1.2273e-16,1.2325e-16,1.2345e-16,1.2345e-16,1.2312e-16,1.2256e-16,1.219e-16,1.2107e-16,1.203e-16,1.1959e-16,1.1907e-16,1.1878e-16,1.1875e-16,1.1895e-16,1.1937e-16,1.2004e-16,1.2083e-16,1.216e-16,1.224e-16,1.2301e-16,1.2345e-16,1.2365e-16,1.2364e-16,1.2336e-16,1.2289e-16,1.2223e-16,1.2148e-16,1.2069e-16,1.1999e-16,1.1938e-16,1.1903e-16,1.1882e-16,1.1885e-16,1.192e-16,1.1965e-16,1.2037e-16,1.2116e-16,1.2194e-16,1.2265e-16,1.2329e-16,1.2377e-16,1.2393e-16,1.2395e-16,1.2369e-16,1.2313e-16,1.225e-16,1.2174e-16,1.2088e-16,1.2017e-16,1.194e-16,1.1901e-16,1.1866e-16,1.1867e-16,1.1894e-16,1.1931e-16,1.2005e-16,1.2076e-16,1.215e-16,1.2243e-16,1.2309e-16,1.2363e-16,1.2421e-16,1.2423e-16,1.2418e-16,1.2385e-16,1.2326e-16,1.2265e-16,1.2185e-16,1.2097e-16,1.2026e-16,1.1953e-16,1.1894e-16,1.1869e-16,1.1854e-16,1.1873e-16,1.1914e-16,1.1955e-16,1.203e-16,1.2113e-16,1.2182e-16,1.227e-16,1.2334e-16,1.2384e-16,1.2427e-16,1.2434e-16,1.2423e-16,1.24e-16,1.2347e-16,1.2278e-16,1.2204e-16,1.2114e-16,1.204e-16,1.1967e-16,1.1904e-16,1.1872e-16,1.185e-16,1.1841e-16,1.1878e-16,1.1916e-16,1.1972e-16,1.2057e-16,1.2138e-16,1.2211e-16,1.2312e-16,1.2365e-16,1.2415e-16,1.2462e-16,1.247e-16,1.2453e-16,1.2433e-16,1.2374e-16,1.2303e-16,1.2235e-16,1.2147e-16,1.205e-16,1.1988e-16,1.1912e-16,1.1851e-16,1.1842e-16,1.1823e-16,1.1827e-16,1.1877e-16,1.1928e-16,1.198e-16,1.2082e-16,1.216e-16,1.2233e-16,1.2324e-16,1.2402e-16,1.2432e-16,1.2483e-16,1.2504e-16,1.2482e-16,1.2455e-16,1.2426e-16,1.2342e-16,1.2269e-16,1.2204e-16,1.2102e-16,1.2024e-16,1.1962e-16,1.1901e-16,1.1833e-16,1.1834e-16,1.1832e-16,1.1827e-16,1.1882e-16,1.1948e-16,1.1995e-16,1.2071e-16,1.2181e-16,1.2244e-16,1.2315e-16,1.2422e-16,1.2464e-16,1.249e-16,1.2538e-16,1.2552e-16,1.2506e-16,1.2477e-16,1.2447e-16,1.2358e-16,1.2274e-16,1.2221e-16,1.2127e-16,1.2011e-16,1.1961e-16,1.1915e-16,1.1841e-16,1.1805e-16,1.1829e-16,1.1825e-16,1.183e-16,1.1901e-16,1.1975e-16,1.2019e-16,1.2102e-16,1.2228e-16,1.2301e-16,1.2353e-16,1.2463e-16,1.2539e-16,1.2554e-16,1.2578e-16,1.2637e-16,1.2616e-16,1.2548e-16,1.2536e-16,1.2503e-16,1.2393e-16,1.2295e-16,1.2255e-16,1.2166e-16,1.2043e-16,1.1972e-16,1.1947e-16,1.1877e-16,1.1801e-16,1.182e-16,1.1846e-16,1.1826e-16,1.1845e-16,1.1955e-16,1.2025e-16,1.2043e-16,1.2147e-16,1.2283e-16,1.2349e-16,1.2386e-16,1.2492e-16,1.2605e-16,1.2618e-16,1.2618e-16,1.2685e-16,1.2721e-16,1.2648e-16,1.2587e-16,1.2588e-16,1.2552e-16,1.2422e-16,1.2317e-16,1.2296e-16,1.2221e-16,1.2071e-16,1.1994e-16,1.1995e-16,1.1943e-16,1.1834e-16,1.1815e-16,1.1883e-16,1.1893e-16,1.1849e-16,1.1898e-16,1.2038e-16,1.2101e-16,1.2102e-16,1.2185e-16,1.2349e-16,1.244e-16,1.2451e-16,1.2519e-16,1.2654e-16,1.2725e-16,1.2697e-16,1.2701e-16,1.2785e-16,1.2823e-16,1.2724e-16,1.2648e-16,1.2668e-16,1.2668e-16,1.2539e-16,1.239e-16,1.2362e-16,1.236e-16,1.2242e-16,1.2079e-16,1.2009e-16,1.2034e-16,1.1989e-16,1.1858e-16,1.1802e-16,1.1869e-16,1.1925e-16,1.1865e-16,1.183e-16,1.194e-16,1.2068e-16,1.21e-16,1.2092e-16,1.2186e-16,1.2373e-16,1.2483e-16,1.2468e-16,1.2517e-16,1.2673e-16,1.2814e-16,1.2817e-16,1.2756e-16,1.2822e-16,1.2946e-16,1.2967e-16,1.284e-16,1.2766e-16,1.2811e-16,1.2844e-16,1.2722e-16,1.2536e-16,1.2468e-16,1.2501e-16,1.2452e-16,1.2261e-16,1.2101e-16,1.2088e-16,1.2134e-16,1.2066e-16,1.1894e-16,1.1807e-16,1.1869e-16,1.1952e-16,1.1909e-16,1.1791e-16,1.1785e-16,1.1933e-16,1.207e-16,1.2066e-16,1.2e-16,1.2065e-16,1.2262e-16,1.2429e-16,1.2444e-16,1.24e-16,1.2479e-16,1.2684e-16,1.2846e-16,1.2842e-16,1.2776e-16,1.2818e-16,1.2982e-16,1.3098e-16,1.3048e-16,1.2925e-16,1.2912e-16,1.3006e-16,1.3067e-16,1.2969e-16,1.2782e-16,1.2698e-16,1.2743e-16,1.2774e-16,1.2661e-16,1.2443e-16,1.2301e-16,1.2317e-16,1.2362e-16,1.2281e-16,1.2076e-16,1.1923e-16,1.1937e-16,1.2024e-16,1.2028e-16,1.1897e-16,1.1743e-16,1.1764e-16,1.191e-16,1.2016e-16,1.199e-16,1.1875e-16,1.1845e-16,1.198e-16,1.2175e-16,1.2278e-16,1.224e-16,1.2166e-16,1.2224e-16,1.242e-16,1.2622e-16,1.2696e-16,1.2633e-16,1.2578e-16,1.2676e-16,1.2888e-16,1.3062e-16,1.3076e-16,1.2984e-16,1.2906e-16,1.2997e-16,1.317e-16,1.3295e-16,1.325e-16,1.309e-16,1.2982e-16,1.3025e-16,1.3158e-16,1.323e-16,1.3131e-16,1.2924e-16,1.2769e-16,1.2776e-16,1.2868e-16,1.2909e-16,1.2787e-16,1.2558e-16,1.2369e-16,1.2348e-16,1.2431e-16,1.248e-16,1.2382e-16,1.2169e-16,1.1961e-16,1.1934e-16,1.2026e-16,1.2124e-16,1.2103e-16,1.1946e-16,1.1748e-16,1.1686e-16,1.1795e-16,1.1958e-16,1.2036e-16,1.1969e-16,1.1811e-16,1.1723e-16,1.1804e-16,1.1998e-16,1.2177e-16,1.2233e-16,1.2151e-16,1.204e-16,1.2062e-16,1.2239e-16,1.2464e-16,1.2621e-16,1.2639e-16,1.255e-16,1.2484e-16,1.2578e-16,1.2783e-16,1.2995e-16,1.3107e-16,1.3078e-16,1.2955e-16,1.2911e-16,1.3025e-16,1.3226e-16,1.3383e-16,1.3435e-16,1.3339e-16,1.318e-16,1.3121e-16,1.3213e-16,1.3371e-16,1.3478e-16,1.3463e-16,1.3295e-16,1.3089e-16,1.3005e-16,1.3074e-16,1.3178e-16,1.3225e-16,1.3149e-16,1.2927e-16,1.2692e-16,1.2608e-16,1.2643e-16,1.2719e-16,1.2742e-16,1.2651e-16,1.2418e-16,1.2175e-16,1.2084e-16,1.214e-16,1.2229e-16,1.2261e-16,1.2185e-16,1.1987e-16,1.1776e-16,1.17e-16,1.179e-16,1.1925e-16,1.2005e-16,1.1994e-16,1.186e-16,1.1682e-16,1.1615e-16,1.1719e-16,1.1916e-16,1.208e-16,1.2155e-16,1.2109e-16,1.1961e-16,1.1879e-16,1.1963e-16,1.2187e-16,1.2417e-16,1.257e-16,1.2612e-16,1.2524e-16,1.2414e-16,1.2414e-16,1.2583e-16,1.2837e-16,1.3054e-16,1.3172e-16,1.3159e-16,1.3041e-16,1.2942e-16,1.2983e-16,1.3164e-16,1.3399e-16,1.3573e-16,1.3635e-16,1.3556e-16,1.3392e-16,1.3267e-16,1.3289e-16,1.3444e-16,1.364e-16,1.3755e-16,1.3747e-16,1.3606e-16,1.3387e-16,1.322e-16,1.3194e-16,1.3309e-16,1.3462e-16,1.3539e-16,1.3501e-16,1.3328e-16,1.3055e-16,1.2828e-16,1.2747e-16,1.2817e-16,1.2943e-16,1.302e-16,1.2986e-16,1.2834e-16,1.2557e-16,1.2287e-16,1.2145e-16,1.2166e-16,1.2291e-16,1.2401e-16,1.2435e-16,1.2337e-16,1.2128e-16,1.1859e-16,1.1662e-16,1.164e-16,1.1761e-16,1.1917e-16,1.2022e-16,1.2047e-16,1.1954e-16,1.1751e-16,1.1541e-16,1.1452e-16,1.1544e-16,1.1735e-16,1.1917e-16,1.2049e-16,1.2095e-16,1.2026e-16,1.1856e-16,1.1717e-16,1.1712e-16,1.1883e-16,1.2119e-16,1.2326e-16,1.2481e-16,1.2569e-16,1.2525e-16,1.2378e-16,1.2268e-16,1.2308e-16,1.25e-16,1.2751e-16,1.2966e-16,1.3133e-16,1.3247e-16,1.323e-16,1.3087e-16,1.2954e-16,1.2964e-16,1.3125e-16,1.3357e-16,1.3545e-16,1.3695e-16,1.3821e-16,1.3836e-16,1.3704e-16,1.3511e-16,1.3414e-16,1.3486e-16,1.3679e-16,1.382e-16,1.3927e-16,1.4025e-16,1.4084e-16,1.4006e-16,1.3788e-16,1.3578e-16,1.35e-16,1.3575e-16,1.3705e-16,1.3778e-16,1.3815e-16,1.3865e-16,1.3862e-16,1.3711e-16,1.3442e-16,1.3213e-16,1.3127e-16,1.3188e-16,1.3278e-16,1.3296e-16,1.3282e-16,1.3298e-16,1.3262e-16,1.3074e-16,1.2773e-16,1.2528e-16,1.2433e-16,1.2487e-16,1.2574e-16,1.2587e-16,1.2553e-16,1.2563e-16,1.2545e-16,1.2381e-16,1.2094e-16,1.1843e-16,1.1734e-16,1.1785e-16,1.1905e-16,1.1965e-16,1.196e-16,1.198e-16,1.2024e-16,1.1958e-16,1.1738e-16,1.1492e-16,1.1361e-16,1.1385e-16,1.1552e-16,1.1687e-16,1.1753e-16,1.179e-16,1.1878e-16,1.1947e-16,1.1878e-16,1.167e-16,1.1512e-16,1.1474e-16,1.159e-16,1.1802e-16,1.2007e-16,1.2097e-16,1.2174e-16,1.2304e-16,1.2399e-16,1.2335e-16,1.2175e-16,1.2064e-16,1.2078e-16,1.2238e-16,1.2498e-16,1.2712e-16,1.2868e-16,1.2955e-16,1.31e-16,1.3207e-16,1.3178e-16,1.3018e-16,1.2903e-16,1.2903e-16,1.3037e-16,1.3271e-16,1.3529e-16,1.3739e-16,1.3813e-16,1.3936e-16,1.4041e-16,1.403e-16,1.388e-16,1.3692e-16,1.3605e-16,1.3645e-16,1.3798e-16,1.4026e-16,1.4217e-16,1.4351e-16,1.4438e-16,1.452e-16,1.4548e-16,1.4439e-16,1.4197e-16,1.3977e-16,1.387e-16,1.3871e-16,1.3974e-16,1.4128e-16,1.4279e-16,1.4399e-16,1.4452e-16,1.4477e-16,1.4441e-16,1.4275e-16,1.399e-16,1.3714e-16,1.355e-16},
      std::array<float, 1001>{ 2.6976e-19,2.6977e-19,2.6978e-19,2.6979e-19,2.6981e-19,2.6982e-19,2.6983e-19,2.6984e-19,2.6985e-19,2.6987e-19,2.6988e-19,2.6989e-19,2.6991e-19,2.6992e-19,2.6994e-19,2.6995e-19,2.6996e-19,2.6997e-19,2.6998e-19,2.6999e-19,2.7001e-19,2.7002e-19,2.7003e-19,2.7004e-19,2.7006e-19,2.7007e-19,2.7009e-19,2.701e-19,2.7012e-19,2.7013e-19,2.7014e-19,2.7015e-19,2.7016e-19,2.7018e-19,2.7019e-19,2.702e-19,2.7021e-19,2.7022e-19,2.7024e-19,2.7025e-19,2.7026e-19,2.7028e-19,2.703e-19,2.7031e-19,2.7033e-19,2.7034e-19,2.7035e-19,2.7036e-19,2.7038e-19,2.7039e-19,2.704e-19,2.7041e-19,2.7042e-19,2.7043e-19,2.7045e-19,2.7046e-19,2.7048e-19,2.7049e-19,2.7051e-19,2.7053e-19,2.7054e-19,2.7056e-19,2.7057e-19,2.7058e-19,2.7059e-19,2.7061e-19,2.7062e-19,2.7063e-19,2.7064e-19,2.7065e-19,2.7067e-19,2.7068e-19,2.707e-19,2.7071e-19,2.7073e-19,2.7075e-19,2.7076e-19,2.7078e-19,2.7079e-19,2.7081e-19,2.7082e-19,2.7083e-19,2.7085e-19,2.7086e-19,2.7087e-19,2.7088e-19,2.7089e-19,2.7091e-19,2.7092e-19,2.7094e-19,2.7095e-19,2.7097e-19,2.7099e-19,2.7101e-19,2.7102e-19,2.7104e-19,2.7105e-19,2.7107e-19,2.7108e-19,2.7109e-19,2.711e-19,2.7112e-19,2.7113e-19,2.7114e-19,2.7115e-19,2.7117e-19,2.7118e-19,2.712e-19,2.7122e-19,2.7123e-19,2.7125e-19,2.7127e-19,2.7129e-19,2.7131e-19,2.7132e-19,2.7134e-19,2.7135e-19,2.7136e-19,2.7137e-19,2.7138e-19,2.7139e-19,2.7141e-19,2.7142e-19,2.7143e-19,2.7145e-19,2.7146e-19,2.7148e-19,2.715e-19,2.7152e-19,2.7154e-19,2.7156e-19,2.7158e-19,2.7159e-19,2.7161e-19,2.7163e-19,2.7164e-19,2.7165e-19,2.7166e-19,2.7168e-19,2.7169e-19,2.717e-19,2.7171e-19,2.7172e-19,2.7174e-19,2.7175e-19,2.7177e-19,2.7179e-19,2.7181e-19,2.7183e-19,2.7185e-19,2.7187e-19,2.7189e-19,2.7191e-19,2.7193e-19,2.7194e-19,2.7196e-19,2.7197e-19,2.7198e-19,2.7199e-19,2.72e-19,2.7201e-19,2.7203e-19,2.7204e-19,2.7206e-19,2.7207e-19,2.7209e-19,2.7211e-19,2.7213e-19,2.7215e-19,2.7218e-19,2.722e-19,2.7222e-19,2.7224e-19,2.7226e-19,2.7227e-19,2.7229e-19,2.723e-19,2.7232e-19,2.7233e-19,2.7234e-19,2.7235e-19,2.7236e-19,2.7237e-19,2.7239e-19,2.724e-19,2.7242e-19,2.7244e-19,2.7246e-19,2.7249e-19,2.7251e-19,2.7253e-19,2.7255e-19,2.7258e-19,2.726e-19,2.7262e-19,2.7264e-19,2.7265e-19,2.7267e-19,2.7268e-19,2.7269e-19,2.727e-19,2.7271e-19,2.7272e-19,2.7274e-19,2.7275e-19,2.7276e-19,2.7278e-19,2.728e-19,2.7282e-19,2.7285e-19,2.7287e-19,2.7289e-19,2.7292e-19,2.7294e-19,2.7296e-19,2.7299e-19,2.7301e-19,2.7303e-19,2.7305e-19,2.7307e-19,2.7308e-19,2.731e-19,2.7311e-19,2.7312e-19,2.7313e-19,2.7314e-19,2.7315e-19,2.7317e-19,2.7318e-19,2.732e-19,2.7322e-19,2.7325e-19,2.7327e-19,2.7329e-19,2.7332e-19,2.7334e-19,2.7336e-19,2.7339e-19,2.7341e-19,2.7344e-19,2.7346e-19,2.7349e-19,2.7351e-19,2.7352e-19,2.7354e-19,2.7355e-19,2.7356e-19,2.7357e-19,2.7358e-19,2.7359e-19,2.7361e-19,2.7362e-19,2.7364e-19,2.7367e-19,2.7369e-19,2.7371e-19,2.7373e-19,2.7376e-19,2.7378e-19,2.738e-19,2.7383e-19,2.7386e-19,2.7388e-19,2.7391e-19,2.7394e-19,2.7396e-19,2.7399e-19,2.7401e-19,2.7402e-19,2.7403e-19,2.7404e-19,2.7405e-19,2.7406e-19,2.7408e-19,2.7409e-19,2.7411e-19,2.7413e-19,2.7415e-19,2.7417e-19,2.742e-19,2.7422e-19,2.7424e-19,2.7426e-19,2.7429e-19,2.7431e-19,2.7434e-19,2.7436e-19,2.7439e-19,2.7442e-19,2.7445e-19,2.7448e-19,2.7451e-19,2.7454e-19,2.7456e-19,2.7458e-19,2.7459e-19,2.746e-19,2.7461e-19,2.7462e-19,2.7464e-19,2.7465e-19,2.7467e-19,2.7469e-19,2.7471e-19,2.7473e-19,2.7476e-19,2.7478e-19,2.748e-19,2.7482e-19,2.7484e-19,2.7486e-19,2.7488e-19,2.7491e-19,2.7493e-19,2.7496e-19,2.7499e-19,2.7503e-19,2.7506e-19,2.7509e-19,2.7512e-19,2.7515e-19,2.7517e-19,2.7519e-19,2.7521e-19,2.7522e-19,2.7523e-19,2.7525e-19,2.7526e-19,2.7528e-19,2.753e-19,2.7532e-19,2.7535e-19,2.7537e-19,2.754e-19,2.7542e-19,2.7544e-19,2.7546e-19,2.7548e-19,2.755e-19,2.7551e-19,2.7553e-19,2.7555e-19,2.7558e-19,2.7561e-19,2.7564e-19,2.7567e-19,2.757e-19,2.7574e-19,2.7577e-19,2.758e-19,2.7583e-19,2.7585e-19,2.7587e-19,2.7589e-19,2.7591e-19,2.7592e-19,2.7594e-19,2.7596e-19,2.7598e-19,2.76e-19,2.7603e-19,2.7606e-19,2.7609e-19,2.7611e-19,2.7614e-19,2.7616e-19,2.7618e-19,2.762e-19,2.7621e-19,2.7623e-19,2.7624e-19,2.7626e-19,2.7628e-19,2.763e-19,2.7633e-19,2.7636e-19,2.7639e-19,2.7642e-19,2.7646e-19,2.7649e-19,2.7652e-19,2.7656e-19,2.7659e-19,2.7661e-19,2.7664e-19,2.7666e-19,2.7668e-19,2.767e-19,2.7672e-19,2.7674e-19,2.7676e-19,2.7679e-19,2.7681e-19,2.7684e-19,2.7687e-19,2.769e-19,2.7693e-19,2.7696e-19,2.7698e-19,2.7699e-19,2.7701e-19,2.7702e-19,2.7703e-19,2.7704e-19,2.7705e-19,2.7707e-19,2.7709e-19,2.7711e-19,2.7714e-19,2.7716e-19,2.7719e-19,2.7723e-19,2.7726e-19,2.7729e-19,2.7733e-19,2.7736e-19,2.7739e-19,2.7743e-19,2.7746e-19,2.7748e-19,2.7751e-19,2.7753e-19,2.7756e-19,2.7758e-19,2.7761e-19,2.7763e-19,2.7766e-19,2.7769e-19,2.7772e-19,2.7775e-19,2.7778e-19,2.7781e-19,2.7783e-19,2.7786e-19,2.7787e-19,2.7788e-19,2.7789e-19,2.779e-19,2.7791e-19,2.7791e-19,2.7792e-19,2.7793e-19,2.7794e-19,2.7796e-19,2.7798e-19,2.78e-19,2.7802e-19,2.7805e-19,2.7808e-19,2.7811e-19,2.7814e-19,2.7817e-19,2.782e-19,2.7823e-19,2.7826e-19,2.7829e-19,2.7832e-19,2.7834e-19,2.7838e-19,2.784e-19,2.7843e-19,2.7845e-19,2.7848e-19,2.7851e-19,2.7853e-19,2.7856e-19,2.7859e-19,2.7862e-19,2.7865e-19,2.7867e-19,2.787e-19,2.7872e-19,2.7874e-19,2.7877e-19,2.7876e-19,2.7876e-19,2.7877e-19,2.7876e-19,2.7877e-19,2.7876e-19,2.7876e-19,2.7876e-19,2.7877e-19,2.7877e-19,2.7878e-19,2.7879e-19,2.788e-19,2.7881e-19,2.7883e-19,2.7884e-19,2.7887e-19,2.7887e-19,2.7889e-19,2.7891e-19,2.7893e-19,2.7894e-19,2.7896e-19,2.7898e-19,2.7899e-19,2.7901e-19,2.7902e-19,2.7904e-19,2.7905e-19,2.7906e-19,2.7908e-19,2.7909e-19,2.791e-19,2.7912e-19,2.7913e-19,2.7914e-19,2.7915e-19,2.7916e-19,2.7916e-19,2.7917e-19,2.7917e-19,2.7916e-19,2.7915e-19,2.7914e-19,2.7912e-19,2.791e-19,2.7908e-19,2.7906e-19,2.7904e-19,2.7901e-19,2.7899e-19,2.7897e-19,2.7895e-19,2.7892e-19,2.789e-19,2.7888e-19,2.7886e-19,2.7884e-19,2.7883e-19,2.7881e-19,2.7879e-19,2.7877e-19,2.7875e-19,2.7873e-19,2.7871e-19,2.7869e-19,2.7867e-19,2.7865e-19,2.7862e-19,2.7859e-19,2.7857e-19,2.7854e-19,2.785e-19,2.7846e-19,2.7842e-19,2.7839e-19,2.7834e-19,2.783e-19,2.7825e-19,2.7821e-19,2.7816e-19,2.781e-19,2.7805e-19,2.7799e-19,2.7793e-19,2.7786e-19,2.7779e-19,2.7772e-19,2.7764e-19,2.7756e-19,2.7747e-19,2.7737e-19,2.7727e-19,2.7717e-19,2.7706e-19,2.7695e-19,2.7684e-19,2.7672e-19,2.766e-19,2.7647e-19,2.7635e-19,2.7622e-19,2.761e-19,2.7596e-19,2.7582e-19,2.7568e-19,2.7555e-19,2.754e-19,2.7526e-19,2.751e-19,2.7495e-19,2.7479e-19,2.7463e-19,2.7447e-19,2.743e-19,2.7413e-19,2.7395e-19,2.7376e-19,2.7357e-19,2.7338e-19,2.7319e-19,2.7299e-19,2.7278e-19,2.7257e-19,2.7235e-19,2.7213e-19,2.7191e-19,2.7167e-19,2.7142e-19,2.7115e-19,2.7087e-19,2.7058e-19,2.7028e-19,2.6998e-19,2.6966e-19,2.6934e-19,2.6901e-19,2.6865e-19,2.6829e-19,2.6793e-19,2.6757e-19,2.6719e-19,2.6678e-19,2.6635e-19,2.6593e-19,2.6551e-19,2.6507e-19,2.6462e-19,2.6412e-19,2.6363e-19,2.6314e-19,2.6265e-19,2.6212e-19,2.6157e-19,2.6101e-19,2.6047e-19,2.5992e-19,2.5932e-19,2.5868e-19,2.5804e-19,2.5743e-19,2.5683e-19,2.5619e-19,2.5548e-19,2.5473e-19,2.5401e-19,2.5333e-19,2.5264e-19,2.5187e-19,2.5104e-19,2.502e-19,2.4941e-19,2.4863e-19,2.4779e-19,2.4687e-19,2.4592e-19,2.4501e-19,2.4414e-19,2.4324e-19,2.4222e-19,2.4171e-19,2.4138e-19,2.4116e-19,2.4099e-19,2.4077e-19,2.4042e-19,2.4002e-19,2.3968e-19,2.3947e-19,2.393e-19,2.3907e-19,2.3871e-19,2.3831e-19,2.3797e-19,2.3775e-19,2.3756e-19,2.3729e-19,2.3691e-19,2.3651e-19,2.362e-19,2.3599e-19,2.3579e-19,2.3547e-19,2.3504e-19,2.3461e-19,2.3429e-19,2.3411e-19,2.3394e-19,2.3363e-19,2.3316e-19,2.3266e-19,2.3229e-19,2.3209e-19,2.3195e-19,2.3168e-19,2.3122e-19,2.3069e-19,2.3026e-19,2.3001e-19,2.2986e-19,2.2962e-19,2.2918e-19,2.2864e-19,2.2817e-19,2.279e-19,2.2773e-19,2.2708e-19,2.2621e-19,2.2518e-19,2.2418e-19,2.234e-19,2.2282e-19,2.2222e-19,2.2137e-19,2.2025e-19,2.1906e-19,2.1808e-19,2.174e-19,2.1685e-19,2.1612e-19,2.1503e-19,2.1373e-19,2.125e-19,2.1158e-19,2.1094e-19,2.1029e-19,2.0933e-19,2.0803e-19,2.0663e-19,2.0545e-19,2.0463e-19,2.0398e-19,2.0313e-19,2.0189e-19,2.0037e-19,1.9891e-19,1.9781e-19,1.9709e-19,1.964e-19,1.9537e-19,1.9389e-19,1.9215e-19,1.9059e-19,1.895e-19,1.8879e-19,1.8805e-19,1.8691e-19,1.8567e-19,1.844e-19,1.8337e-19,1.8285e-19,1.8268e-19,1.8248e-19,1.8194e-19,1.809e-19,1.7964e-19,1.7863e-19,1.7808e-19,1.7783e-19,1.7759e-19,1.7702e-19,1.7594e-19,1.7465e-19,1.7361e-19,1.7301e-19,1.7275e-19,1.7257e-19,1.7206e-19,1.7099e-19,1.6962e-19,1.6841e-19,1.6761e-19,1.6725e-19,1.6716e-19,1.6685e-19,1.6597e-19,1.6465e-19,1.6327e-19,1.6215e-19,1.6151e-19,1.6132e-19,1.6117e-19,1.6058e-19,1.5948e-19,1.5811e-19,1.5728e-19,1.5691e-19,1.571e-19,1.5757e-19,1.5783e-19,1.5762e-19,1.5698e-19,1.5615e-19,1.555e-19,1.5542e-19,1.5582e-19,1.563e-19,1.5644e-19,1.5606e-19,1.5521e-19,1.5426e-19,1.5371e-19,1.5381e-19,1.5434e-19,1.5491e-19,1.5506e-19,1.5459e-19,1.536e-19,1.5257e-19,1.5203e-19,1.5212e-19,1.527e-19,1.5333e-19,1.5353e-19,1.5305e-19,1.5207e-19,1.5106e-19,1.5048e-19,1.5053e-19,1.5106e-19,1.5135e-19,1.5121e-19,1.5042e-19,1.4915e-19,1.4781e-19,1.4685e-19,1.4654e-19,1.468e-19,1.4713e-19,1.4704e-19,1.4631e-19,1.4505e-19,1.436e-19,1.4242e-19,1.4196e-19,1.422e-19,1.4262e-19,1.4273e-19,1.4228e-19,1.4122e-19,1.397e-19,1.3821e-19,1.3739e-19,1.3737e-19,1.3769e-19,1.3794e-19,1.3786e-19,1.3726e-19,1.3598e-19,1.3436e-19,1.3316e-19,1.3276e-19,1.3316e-19,1.3369e-19,1.3416e-19,1.3439e-19,1.3398e-19,1.3291e-19,1.3186e-19,1.3148e-19,1.3166e-19,1.3199e-19,1.3243e-19,1.3293e-19,1.3309e-19,1.3243e-19,1.3124e-19,1.304e-19,1.3021e-19,1.3029e-19,1.3051e-19,1.3099e-19,1.3167e-19,1.3187e-19,1.3116e-19,1.3001e-19,1.2933e-19,1.2906e-19,1.2883e-19,1.2881e-19,1.2925e-19,1.3008e-19,1.2845e-19,1.2571e-19,1.2254e-19,1.1976e-19,1.1739e-19,1.1465e-19,1.1163e-19,1.0946e-19,1.1071e-19,1.0597e-19,1.0319e-19,9.9782e-20,9.6941e-20,9.4179e-20,9.1042e-20,8.7408e-20,8.4353e-20,8.2344e-20,8.0216e-20,7.7287e-20,7.3503e-20,7.016e-20,6.7448e-20,6.4137e-20,5.9912e-20,5.581e-20,5.2759e-20,5.0103e-20,5.0205e-20,4.9996e-20,4.9499e-20,4.9806e-20,5.0005e-20,4.94e-20,4.8058e-20,4.7676e-20,4.7969e-20,4.8285e-20,4.8451e-20,4.8888e-20,4.8514e-20,4.8909e-20,4.8713e-20,4.7814e-20,4.6831e-20,4.657e-20,4.6553e-20,4.65e-20,4.663e-20,4.6898e-20,4.7319e-20,4.7737e-20,4.7627e-20,4.6429e-20,4.5798e-20,4.5533e-20,4.5221e-20,4.4861e-20,4.5659e-20,4.5407e-20,4.6039e-20,4.6524e-20,4.6233e-20,4.5356e-20,4.4913e-20,4.4514e-20,4.411e-20,4.3363e-20,4.3273e-20,4.3597e-20,4.4833e-20,4.5304e-20,4.5212e-20,4.456e-20,4.399e-20,4.3704e-20,4.3205e-20,4.2218e-20,4.1518e-20,4.1596e-20,4.2525e-20,4.3309e-20,4.3769e-20,4.3278e-20,4.2747e-20,4.2501e-20,4.2339e-20,4.1291e-20,3.993e-20,3.9266e-20,3.9513e-20,4.0456e-20,4.1335e-20,4.1432e-20,4.0809e-20,4.0606e-20,4.075e-20,4.0612e-20,3.9373e-20,3.7795e-20,3.7154e-20,3.7383e-20,3.8561e-20,3.935e-20,3.9051e-20,3.8475e-20,3.8366e-20,3.8759e-20,3.8919e-20,3.7756e-20,3.5922e-20,3.4938e-20,3.5055e-20,3.6193e-20,3.7232e-20,3.7058e-20,3.625e-20,3.6043e-20,3.6437e-20,3.714e-20,3.6742e-20,3.4879e-20,3.3093e-20,3.2634e-20,3.34e-20} };

    const static float base_density = 2.3864e+16f;//km^-3

    const static float Hp = 30000e-3f; // Continental value

    const static float pb = 2e6f; // background value in km^-3

    template <typename Float, typename UInt32, typename Mask, typename Spectrum, typename Wavelength>
    class PolarAntarticAerosol : public GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength> {
        const float pb_bd;

    public:
        PolarAntarticAerosol() : GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength>(tabulatedValues),
                                 pb_bd(pb / base_density) {}

        Float get_density(const Float &z) const override {
            const Mask msk = z >= Float(0) && z <= Float(86);
            const Float height = enoki::select(
                    msk,
                    z,
                    Float(0.f)
            );
            return enoki::select(
                    msk,
                    Float(base_density) * (enoki::exp(-height / Float(Hp)) + Float(pb_bd)),
                    height
            );
        }

        [[nodiscard]] float get_density_float(const float &z) const override {
            return base_density * (enoki::exp(-z / Hp) + pb_bd);
        }
    };
}
#endif //_AEROSOL_POLAR_ANTARTIC