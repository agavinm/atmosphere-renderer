#ifndef _AEROSOL_RURAL
#define _AEROSOL_RURAL

#include <array>
#include <mitsuba/atmosphere/GlobalAtmosphericAerosol.h>

namespace ruralAerosol {

	// row[0] = wavelength
	// row[1] = cross section absorption coefficient
	// row[2] = cross section scattering coefficient
    const static std::array<std::array<float, 1001>, 3> tabulatedValues =
	{ std::array<float, 1001>{ 100,100.23,100.46,100.69,100.93,101.16,101.39,101.62,101.86,102.09,102.33,102.57,102.8,103.04,103.28,103.51,103.75,103.99,104.23,104.47,104.71,104.95,105.2,105.44,105.68,105.93,106.17,106.41,106.66,106.91,107.15,107.4,107.65,107.89,108.14,108.39,108.64,108.89,109.14,109.4,109.65,109.9,110.15,110.41,110.66,110.92,111.17,111.43,111.69,111.94,112.2,112.46,112.72,112.98,113.24,113.5,113.76,114.03,114.29,114.55,114.82,115.08,115.35,115.61,115.88,116.14,116.41,116.68,116.95,117.22,117.49,117.76,118.03,118.3,118.58,118.85,119.12,119.4,119.67,119.95,120.23,120.5,120.78,121.06,121.34,121.62,121.9,122.18,122.46,122.74,123.03,123.31,123.59,123.88,124.17,124.45,124.74,125.03,125.31,125.6,125.89,126.18,126.47,126.77,127.06,127.35,127.64,127.94,128.23,128.53,128.82,129.12,129.42,129.72,130.02,130.32,130.62,130.92,131.22,131.52,131.83,132.13,132.43,132.74,133.05,133.35,133.66,133.97,134.28,134.59,134.9,135.21,135.52,135.83,136.14,136.46,136.77,137.09,137.4,137.72,138.04,138.36,138.68,139,139.32,139.64,139.96,140.28,140.6,140.93,141.25,141.58,141.91,142.23,142.56,142.89,143.22,143.55,143.88,144.21,144.54,144.88,145.21,145.55,145.88,146.22,146.55,146.89,147.23,147.57,147.91,148.25,148.59,148.94,149.28,149.62,149.97,150.31,150.66,151.01,151.36,151.71,152.05,152.41,152.76,153.11,153.46,153.82,154.17,154.53,154.88,155.24,155.6,155.96,156.31,156.68,157.04,157.4,157.76,158.12,158.49,158.85,159.22,159.59,159.96,160.32,160.69,161.06,161.44,161.81,162.18,162.55,162.93,163.31,163.68,164.06,164.44,164.82,165.2,165.58,165.96,166.34,166.72,167.11,167.49,167.88,168.27,168.66,169.04,169.43,169.82,170.22,170.61,171,171.4,171.79,172.19,172.58,172.98,173.38,173.78,174.18,174.58,174.98,175.39,175.79,176.2,176.6,177.01,177.42,177.83,178.24,178.65,179.06,179.47,179.89,180.3,180.72,181.13,181.55,181.97,182.39,182.81,183.23,183.65,184.08,184.5,184.93,185.35,185.78,186.21,186.64,187.07,187.5,187.93,188.36,188.8,189.23,189.67,190.11,190.55,190.99,191.43,191.87,192.31,192.75,193.2,193.64,194.09,194.54,194.98,195.43,195.88,196.34,196.79,197.24,197.7,198.15,198.61,199.07,199.53,199.99,200.45,200.91,201.37,201.84,202.3,202.77,203.24,203.7,204.17,204.64,205.12,205.59,206.06,206.54,207.01,207.49,207.97,208.45,208.93,209.41,209.89,210.38,210.86,211.35,211.84,212.32,212.81,213.3,213.8,214.29,214.78,215.28,215.77,216.27,216.77,217.27,217.77,218.27,218.78,219.28,219.79,220.29,220.8,221.31,221.82,222.33,222.84,223.36,223.87,224.39,224.91,225.42,225.94,226.46,226.99,227.51,228.03,228.56,229.09,229.61,230.14,230.67,231.21,231.74,232.27,232.81,233.35,233.88,234.42,234.96,235.5,236.05,236.59,237.14,237.68,238.23,238.78,239.33,239.88,240.44,240.99,241.55,242.1,242.66,243.22,243.78,244.34,244.91,245.47,246.04,246.6,247.17,247.74,248.31,248.89,249.46,250.03,250.61,251.19,251.77,252.35,252.93,253.51,254.1,254.68,255.27,255.86,256.45,257.04,257.63,258.23,258.82,259.42,260.02,260.62,261.22,261.82,262.42,263.03,263.63,264.24,264.85,265.46,266.07,266.69,267.3,267.92,268.53,269.15,269.77,270.4,271.02,271.64,272.27,272.9,273.53,274.16,274.79,275.42,276.06,276.69,277.33,277.97,278.61,279.25,279.9,280.54,281.19,281.84,282.49,283.14,283.79,284.45,285.1,285.76,286.42,287.08,287.74,288.4,289.07,289.73,290.4,291.07,291.74,292.42,293.09,293.76,294.44,295.12,295.8,296.48,297.17,297.85,298.54,299.23,299.92,300.61,301.3,302,302.69,303.39,304.09,304.79,305.49,306.2,306.9,307.61,308.32,309.03,309.74,310.46,311.17,311.89,312.61,313.33,314.05,314.77,315.5,316.23,316.96,317.69,318.42,319.15,319.89,320.63,321.37,322.11,322.85,323.59,324.34,325.09,325.84,326.59,327.34,328.1,328.85,329.61,330.37,331.13,331.89,332.66,333.43,334.19,334.97,335.74,336.51,337.29,338.06,338.84,339.63,340.41,341.19,341.98,342.77,343.56,344.35,345.14,345.94,346.74,347.54,348.34,349.14,349.95,350.75,351.56,352.37,353.18,354,354.81,355.63,356.45,357.27,358.1,358.92,359.75,360.58,361.41,362.24,363.08,363.92,364.75,365.59,366.44,367.28,368.13,368.98,369.83,370.68,371.54,372.39,373.25,374.11,374.97,375.84,376.7,377.57,378.44,379.32,380.19,381.07,381.94,382.82,383.71,384.59,385.48,386.37,387.26,388.15,389.05,389.94,390.84,391.74,392.64,393.55,394.46,395.37,396.28,397.19,398.11,399.02,399.94,400.87,401.79,402.72,403.65,404.58,405.51,406.44,407.38,408.32,409.26,410.2,411.15,412.1,413.05,414,414.95,415.91,416.87,417.83,418.79,419.76,420.73,421.7,422.67,423.64,424.62,425.6,426.58,427.56,428.55,429.54,430.53,431.52,432.51,433.51,434.51,435.51,436.52,437.52,438.53,439.54,440.55,441.57,442.59,443.61,444.63,445.66,446.68,447.71,448.75,449.78,450.82,451.86,452.9,453.94,454.99,456.04,457.09,458.14,459.2,460.26,461.32,462.38,463.45,464.52,465.59,466.66,467.74,468.81,469.89,470.98,472.06,473.15,474.24,475.34,476.43,477.53,478.63,479.73,480.84,481.95,483.06,484.17,485.29,486.41,487.53,488.65,489.78,490.91,492.04,493.17,494.31,495.45,496.59,497.74,498.88,500.03,501.19,502.34,503.5,504.66,505.82,506.99,508.16,509.33,510.51,511.68,512.86,514.04,515.23,516.42,517.61,518.8,520,521.19,522.4,523.6,524.81,526.02,527.23,528.45,529.66,530.88,532.11,533.33,534.56,535.8,537.03,538.27,539.51,540.75,542,543.25,544.5,545.76,547.02,548.28,549.54,550.81,552.08,553.35,554.63,555.9,557.19,558.47,559.76,561.05,562.34,563.64,564.94,566.24,567.54,568.85,570.16,571.48,572.8,574.12,575.44,576.77,578.1,579.43,580.76,582.1,583.45,584.79,586.14,587.49,588.84,590.2,591.56,592.93,594.29,595.66,597.04,598.41,599.79,601.17,602.56,603.95,605.34,606.74,608.13,609.54,610.94,612.35,613.76,615.18,616.6,618.02,619.44,620.87,622.3,623.73,625.17,626.61,628.06,629.51,630.96,632.41,633.87,635.33,636.8,638.26,639.73,641.21,642.69,644.17,645.65,647.14,648.63,650.13,651.63,653.13,654.64,656.15,657.66,659.17,660.69,662.22,663.74,665.27,666.81,668.34,669.88,671.43,672.98,674.53,676.08,677.64,679.2,680.77,682.34,683.91,685.49,687.07,688.65,690.24,691.83,693.43,695.02,696.63,698.23,699.84,701.46,703.07,704.69,706.32,707.95,709.58,711.21,712.85,714.5,716.14,717.79,719.45,721.11,722.77,724.44,726.11,727.78,729.46,731.14,732.82,734.51,736.21,737.9,739.61,741.31,743.02,744.73,746.45,748.17,749.89,751.62,753.36,755.09,756.83,758.58,760.33,762.08,763.84,765.6,767.36,769.13,770.9,772.68,774.46,776.25,778.04,779.83,781.63,783.43,785.24,787.05,788.86,790.68,792.5,794.33,796.16,797.99,799.83,801.68,803.53,805.38,807.24,809.1,810.96,812.83,814.7,816.58,818.46,820.35,822.24,824.14,826.04,827.94,829.85,831.76,833.68,835.6,837.53,839.46,841.4,843.33,845.28,847.23,849.18,851.14,853.1,855.07,857.04,859.01,860.99,862.98,864.97,866.96,868.96,870.96,872.97,874.98,877,879.02,881.05,883.08,885.12,887.16,889.2,891.25,893.31,895.36,897.43,899.5,901.57,903.65,905.73,907.82,909.91,912.01,914.11,916.22,918.33,920.45,922.57,924.7,926.83,928.97,931.11,933.25,935.41,937.56,939.72,941.89,944.06,946.24,948.42,950.6,952.8,954.99,957.19,959.4,961.61,963.83,966.05,968.28,970.51,972.75,974.99,977.24,979.49,981.75,984.01,986.28,988.55,990.83,993.12,995.41,997.7,1000},
      std::array<float, 1001>{ 1.6719e-20,1.6653e-20,1.6586e-20,1.6519e-20,1.6451e-20,1.6383e-20,1.6314e-20,1.6245e-20,1.6175e-20,1.6105e-20,1.6035e-20,1.5965e-20,1.5894e-20,1.5823e-20,1.5753e-20,1.5682e-20,1.5612e-20,1.5541e-20,1.5471e-20,1.5401e-20,1.5331e-20,1.5261e-20,1.5192e-20,1.5123e-20,1.5055e-20,1.4987e-20,1.4919e-20,1.4852e-20,1.4786e-20,1.472e-20,1.4654e-20,1.459e-20,1.4526e-20,1.4462e-20,1.4399e-20,1.4337e-20,1.4276e-20,1.4215e-20,1.4156e-20,1.4096e-20,1.4038e-20,1.398e-20,1.3923e-20,1.3867e-20,1.3812e-20,1.3757e-20,1.3703e-20,1.365e-20,1.3598e-20,1.3546e-20,1.3495e-20,1.3445e-20,1.3396e-20,1.3347e-20,1.3299e-20,1.3251e-20,1.3204e-20,1.3158e-20,1.3113e-20,1.3068e-20,1.3024e-20,1.298e-20,1.2937e-20,1.2894e-20,1.2852e-20,1.2811e-20,1.277e-20,1.2729e-20,1.2689e-20,1.265e-20,1.2611e-20,1.2572e-20,1.2534e-20,1.2496e-20,1.2458e-20,1.2421e-20,1.2384e-20,1.2347e-20,1.2311e-20,1.2275e-20,1.2239e-20,1.2203e-20,1.2168e-20,1.2133e-20,1.2098e-20,1.2063e-20,1.2028e-20,1.1993e-20,1.1958e-20,1.1924e-20,1.1889e-20,1.1855e-20,1.182e-20,1.1786e-20,1.1751e-20,1.1717e-20,1.1682e-20,1.1648e-20,1.1613e-20,1.1578e-20,1.1543e-20,1.1508e-20,1.1473e-20,1.1437e-20,1.1402e-20,1.1366e-20,1.133e-20,1.1293e-20,1.1257e-20,1.122e-20,1.1183e-20,1.1146e-20,1.1108e-20,1.1071e-20,1.1032e-20,1.0994e-20,1.0955e-20,1.0916e-20,1.0877e-20,1.0837e-20,1.0797e-20,1.0756e-20,1.0715e-20,1.0674e-20,1.0632e-20,1.059e-20,1.0548e-20,1.0505e-20,1.0462e-20,1.0418e-20,1.0374e-20,1.033e-20,1.0285e-20,1.024e-20,1.0195e-20,1.0149e-20,1.0103e-20,1.0056e-20,1.0009e-20,9.9622e-21,9.9146e-21,9.8666e-21,9.8183e-21,9.7697e-21,9.7207e-21,9.6715e-21,9.6219e-21,9.5721e-21,9.522e-21,9.4716e-21,9.4209e-21,9.3701e-21,9.3189e-21,9.2676e-21,9.2161e-21,9.1643e-21,9.1124e-21,9.0603e-21,9.0081e-21,8.9557e-21,8.9032e-21,8.8506e-21,8.7979e-21,8.7452e-21,8.6923e-21,8.6394e-21,8.5865e-21,8.5335e-21,8.4806e-21,8.4276e-21,8.3747e-21,8.3218e-21,8.269e-21,8.2162e-21,8.1635e-21,8.1109e-21,8.0583e-21,8.0057e-21,7.9532e-21,7.9008e-21,7.8486e-21,7.7965e-21,7.7447e-21,7.693e-21,7.6415e-21,7.5903e-21,7.5392e-21,7.4884e-21,7.4379e-21,7.3876e-21,7.3376e-21,7.2878e-21,7.2383e-21,7.1891e-21,7.1401e-21,7.0915e-21,7.0432e-21,6.9951e-21,6.9474e-21,6.9e-21,6.8529e-21,6.8062e-21,6.7597e-21,6.7136e-21,6.6678e-21,6.6224e-21,6.5773e-21,6.5326e-21,6.4881e-21,6.4441e-21,6.4003e-21,6.3569e-21,6.3139e-21,6.2712e-21,6.2288e-21,6.1868e-21,6.1452e-21,6.1038e-21,6.0629e-21,6.0222e-21,5.9819e-21,5.9419e-21,5.9023e-21,5.863e-21,5.824e-21,5.7854e-21,5.7471e-21,5.7091e-21,5.6714e-21,5.6341e-21,5.5971e-21,5.5603e-21,5.5239e-21,5.4878e-21,5.452e-21,5.4165e-21,5.3813e-21,5.3464e-21,5.3117e-21,5.2774e-21,5.2433e-21,5.2096e-21,5.1761e-21,5.1428e-21,5.1098e-21,5.0771e-21,5.0447e-21,5.0125e-21,4.9806e-21,4.9489e-21,4.9174e-21,4.8862e-21,4.8553e-21,4.8245e-21,4.794e-21,4.7638e-21,4.7337e-21,4.7039e-21,4.6743e-21,4.6449e-21,4.6157e-21,4.5867e-21,4.558e-21,4.5294e-21,4.501e-21,4.4728e-21,4.4448e-21,4.417e-21,4.3894e-21,4.362e-21,4.3347e-21,4.3077e-21,4.2807e-21,4.254e-21,4.2274e-21,4.201e-21,4.1748e-21,4.1487e-21,4.1228e-21,4.0971e-21,4.0714e-21,4.046e-21,4.0207e-21,3.9955e-21,3.9705e-21,3.9456e-21,3.9209e-21,3.8963e-21,3.8718e-21,3.8475e-21,3.8233e-21,3.7993e-21,3.7753e-21,3.7515e-21,3.7279e-21,3.7043e-21,3.6809e-21,3.6576e-21,3.6344e-21,3.6113e-21,3.5884e-21,3.5656e-21,3.5428e-21,3.5201e-21,3.4976e-21,3.4751e-21,3.4528e-21,3.4306e-21,3.4085e-21,3.3865e-21,3.3646e-21,3.3428e-21,3.3211e-21,3.2996e-21,3.2781e-21,3.2567e-21,3.2355e-21,3.2143e-21,3.1932e-21,3.1723e-21,3.1514e-21,3.1307e-21,3.11e-21,3.0895e-21,3.069e-21,3.0486e-21,3.0284e-21,3.0082e-21,2.9881e-21,2.9681e-21,2.9482e-21,2.9284e-21,2.9087e-21,2.8891e-21,2.8696e-21,2.8502e-21,2.8309e-21,2.8116e-21,2.7925e-21,2.7734e-21,2.7545e-21,2.7356e-21,2.7168e-21,2.6981e-21,2.6795e-21,2.661e-21,2.6426e-21,2.6242e-21,2.606e-21,2.5878e-21,2.5698e-21,2.5518e-21,2.5339e-21,2.5161e-21,2.4984e-21,2.4808e-21,2.4633e-21,2.4458e-21,2.4284e-21,2.4112e-21,2.394e-21,2.3769e-21,2.3599e-21,2.343e-21,2.3261e-21,2.3094e-21,2.2927e-21,2.2761e-21,2.2596e-21,2.2432e-21,2.2269e-21,2.2107e-21,2.1946e-21,2.1785e-21,2.1625e-21,2.1466e-21,2.1308e-21,2.1151e-21,2.0995e-21,2.0839e-21,2.0685e-21,2.0531e-21,2.0378e-21,2.0226e-21,2.0075e-21,1.9925e-21,1.9775e-21,1.9626e-21,1.9479e-21,1.9332e-21,1.9186e-21,1.904e-21,1.8896e-21,1.8752e-21,1.8609e-21,1.8467e-21,1.8326e-21,1.8186e-21,1.8047e-21,1.7908e-21,1.777e-21,1.7632e-21,1.7496e-21,1.736e-21,1.7226e-21,1.7092e-21,1.6959e-21,1.6826e-21,1.6695e-21,1.6564e-21,1.6434e-21,1.6305e-21,1.6177e-21,1.605e-21,1.5923e-21,1.5797e-21,1.5672e-21,1.5548e-21,1.5425e-21,1.5302e-21,1.518e-21,1.5059e-21,1.4939e-21,1.482e-21,1.4701e-21,1.4583e-21,1.4466e-21,1.435e-21,1.4234e-21,1.412e-21,1.4006e-21,1.3892e-21,1.378e-21,1.3668e-21,1.3557e-21,1.3447e-21,1.3337e-21,1.3229e-21,1.3121e-21,1.3013e-21,1.2907e-21,1.2801e-21,1.2696e-21,1.2592e-21,1.2488e-21,1.2386e-21,1.2283e-21,1.2182e-21,1.2081e-21,1.1981e-21,1.1882e-21,1.1784e-21,1.1686e-21,1.1589e-21,1.1492e-21,1.1397e-21,1.1301e-21,1.1207e-21,1.1113e-21,1.102e-21,1.0928e-21,1.0837e-21,1.0746e-21,1.0655e-21,1.0566e-21,1.0477e-21,1.0388e-21,1.0301e-21,1.0214e-21,1.0127e-21,1.0042e-21,9.9566e-22,9.8722e-22,9.7884e-22,9.7052e-22,9.6227e-22,9.5408e-22,9.4595e-22,9.3788e-22,9.2987e-22,9.2192e-22,9.1403e-22,9.062e-22,8.9843e-22,8.9072e-22,8.8307e-22,8.7548e-22,8.6795e-22,8.6048e-22,8.5306e-22,8.4571e-22,8.3841e-22,8.3117e-22,8.2399e-22,8.1686e-22,8.0979e-22,8.0277e-22,7.9581e-22,7.8891e-22,7.8206e-22,7.7527e-22,7.6852e-22,7.6184e-22,7.5521e-22,7.4863e-22,7.421e-22,7.3562e-22,7.292e-22,7.2283e-22,7.1651e-22,7.1024e-22,7.0403e-22,6.9786e-22,6.9174e-22,6.8568e-22,6.7966e-22,6.7369e-22,6.6777e-22,6.619e-22,6.5608e-22,6.5031e-22,6.4458e-22,6.389e-22,6.3327e-22,6.2768e-22,6.2214e-22,6.1664e-22,6.112e-22,6.0579e-22,6.0043e-22,5.9512e-22,5.8985e-22,5.8462e-22,5.7944e-22,5.743e-22,5.6921e-22,5.6415e-22,5.5914e-22,5.5417e-22,5.4925e-22,5.4436e-22,5.3952e-22,5.3471e-22,5.2995e-22,5.2523e-22,5.2054e-22,5.1589e-22,5.1129e-22,5.0672e-22,5.0219e-22,4.977e-22,4.9325e-22,4.8883e-22,4.8446e-22,4.8012e-22,4.7582e-22,4.7155e-22,4.6732e-22,4.6313e-22,4.5898e-22,4.5486e-22,4.5077e-22,4.4672e-22,4.4271e-22,4.3873e-22,4.3478e-22,4.3087e-22,4.2699e-22,4.2315e-22,4.1934e-22,4.1556e-22,4.1182e-22,4.081e-22,4.0442e-22,4.0078e-22,3.9716e-22,3.9358e-22,3.9002e-22,3.865e-22,3.8301e-22,3.7955e-22,3.7612e-22,3.7272e-22,3.6935e-22,3.66e-22,3.6269e-22,3.5941e-22,3.5615e-22,3.5293e-22,3.4973e-22,3.4656e-22,3.4342e-22,3.4031e-22,3.3722e-22,3.3416e-22,3.3113e-22,3.2813e-22,3.2515e-22,3.222e-22,3.1927e-22,3.1637e-22,3.1349e-22,3.1065e-22,3.0782e-22,3.0502e-22,3.0224e-22,2.9949e-22,2.9677e-22,2.9406e-22,2.9139e-22,2.8873e-22,2.861e-22,2.835e-22,2.8091e-22,2.7835e-22,2.7581e-22,2.733e-22,2.7081e-22,2.6833e-22,2.6589e-22,2.6346e-22,2.6106e-22,2.5867e-22,2.5631e-22,2.5397e-22,2.5165e-22,2.4935e-22,2.4707e-22,2.4482e-22,2.4258e-22,2.4036e-22,2.3816e-22,2.3599e-22,2.3383e-22,2.3169e-22,2.2957e-22,2.2747e-22,2.2539e-22,2.2333e-22,2.2128e-22,2.1926e-22,2.1725e-22,2.1526e-22,2.1329e-22,2.1133e-22,2.094e-22,2.0748e-22,2.0558e-22,2.037e-22,2.0183e-22,1.9998e-22,1.9815e-22,1.9633e-22,1.9453e-22,1.9274e-22,1.9097e-22,1.8922e-22,1.8748e-22,1.8576e-22,1.8405e-22,1.8236e-22,1.8069e-22,1.7903e-22,1.7738e-22,1.7575e-22,1.7414e-22,1.7254e-22,1.7095e-22,1.6938e-22,1.6783e-22,1.6628e-22,1.6476e-22,1.6324e-22,1.6174e-22,1.6025e-22,1.5878e-22,1.5732e-22,1.5587e-22,1.5444e-22,1.5302e-22,1.5161e-22,1.5022e-22,1.4884e-22,1.4747e-22,1.4611e-22,1.4477e-22,1.4344e-22,1.4212e-22,1.4081e-22,1.3952e-22,1.3823e-22,1.3696e-22,1.357e-22,1.3445e-22,1.3322e-22,1.3199e-22,1.3078e-22,1.2958e-22,1.2838e-22,1.272e-22,1.2603e-22,1.2487e-22,1.2372e-22,1.2258e-22,1.2145e-22,1.2033e-22,1.1922e-22,1.1812e-22,1.1703e-22,1.1595e-22,1.1488e-22,1.1382e-22,1.1277e-22,1.1173e-22,1.107e-22,1.0968e-22,1.0866e-22,1.0766e-22,1.0667e-22,1.0568e-22,1.0471e-22,1.0374e-22,1.0278e-22,1.0183e-22,1.0089e-22,9.9963e-23,9.904e-23,9.8126e-23,9.722e-23,9.6322e-23,9.5433e-23,9.4551e-23,9.3678e-23,9.2813e-23,9.1956e-23,9.1107e-23,9.0265e-23,8.9431e-23,8.8605e-23,8.7787e-23,8.6976e-23,8.6173e-23,8.5379e-23,8.4592e-23,8.3813e-23,8.304e-23,8.2275e-23,8.1516e-23,8.0765e-23,8.0021e-23,7.9283e-23,7.8553e-23,7.7829e-23,7.7111e-23,7.6401e-23,7.5697e-23,7.4999e-23,7.4308e-23,7.3623e-23,7.2944e-23,7.2272e-23,7.1606e-23,7.0946e-23,7.0292e-23,6.9645e-23,6.9003e-23,6.8367e-23,6.7737e-23,6.7113e-23,6.6494e-23,6.5881e-23,6.5274e-23,6.4673e-23,6.4077e-23,6.3486e-23,6.2901e-23,6.2322e-23,6.1748e-23,6.1179e-23,6.0616e-23,6.0059e-23,5.9506e-23,5.8959e-23,5.8417e-23,5.788e-23,5.7348e-23,5.6821e-23,5.6298e-23,5.5781e-23,5.5268e-23,5.476e-23,5.4256e-23,5.3758e-23,5.3264e-23,5.2774e-23,5.2289e-23,5.1808e-23,5.1332e-23,5.086e-23,5.0393e-23,4.9929e-23,4.9471e-23,4.9016e-23,4.8565e-23,4.8119e-23,4.7677e-23,4.7239e-23,4.6804e-23,4.6374e-23,4.5948e-23,4.5526e-23,4.5108e-23,4.4693e-23,4.4282e-23,4.3874e-23,4.347e-23,4.3069e-23,4.2672e-23,4.2279e-23,4.189e-23,4.1504e-23,4.1121e-23,4.0742e-23,4.0367e-23,3.9995e-23,3.9627e-23,3.9261e-23,3.89e-23,3.8541e-23,3.8186e-23,3.7834e-23,3.7486e-23,3.714e-23,3.6798e-23,3.6459e-23,3.6123e-23,3.5791e-23,3.5461e-23,3.5134e-23,3.4811e-23,3.449e-23,3.4172e-23,3.3857e-23,3.3546e-23,3.3237e-23,3.293e-23,3.2628e-23,3.2329e-23,3.2032e-23,3.1739e-23,3.1448e-23,3.1159e-23,3.0873e-23,3.059e-23,3.031e-23,3.0032e-23,2.9757e-23,2.9484e-23,2.9214e-23,2.8946e-23,2.8681e-23,2.8418e-23,2.8157e-23,2.7899e-23,2.7644e-23,2.739e-23,2.7139e-23,2.6891e-23,2.6644e-23,2.64e-23,2.6158e-23,2.5919e-23,2.5681e-23,2.5446e-23,2.5213e-23,2.4982e-23,2.4735e-23,2.449e-23,2.4248e-23,2.4008e-23,2.3771e-23,2.3536e-23,2.3304e-23,2.3074e-23,2.2847e-23,2.2622e-23,2.24e-23,2.218e-23,2.1962e-23,2.1747e-23,2.1534e-23,2.1324e-23,2.1116e-23,2.091e-23,2.0706e-23,2.0504e-23,2.0305e-23,2.0108e-23,1.9913e-23,1.972e-23,1.9529e-23,1.934e-23,1.9153e-23,1.8968e-23,1.8793e-23,1.862e-23,1.8448e-23,1.8278e-23,1.8109e-23,1.7942e-23,1.7777e-23,1.7613e-23,1.7451e-23,1.729e-23,1.713e-23,1.6972e-23,1.6816e-23,1.6661e-23,1.6507e-23,1.6355e-23,1.6204e-23,1.6055e-23,1.5907e-23,1.576e-23,1.5615e-23,1.5471e-23,1.5328e-23,1.5187e-23,1.5047e-23,1.4908e-23,1.4771e-23,1.4635e-23,1.45e-23,1.4366e-23,1.4234e-23,1.4102e-23,1.3972e-23,1.3844e-23,1.3716e-23,1.359e-23,1.3464e-23,1.334e-23,1.3217e-23,1.3095e-23,1.2975e-23,1.2855e-23,1.2736e-23,1.2619e-23,1.2503e-23,1.2387e-23,1.2273e-23,1.216e-23,1.2048e-23,1.1937e-23,1.1827e-23,1.1734e-23,1.1646e-23,1.1559e-23,1.1473e-23,1.1388e-23,1.1303e-23,1.1219e-23,1.1136e-23,1.1053e-23,1.0971e-23,1.089e-23,1.081e-23,1.073e-23,1.065e-23,1.0572e-23,1.0494e-23,1.0416e-23,1.034e-23,1.0264e-23,1.0188e-23,1.0113e-23,1.0039e-23,9.965e-24,9.8918e-24,9.8192e-24,9.7472e-24,9.6757e-24,9.6048e-24,9.5345e-24,9.4647e-24,9.3954e-24,9.3267e-24,9.2585e-24,9.1908e-24,9.1236e-24,9.057e-24,8.9908e-24,8.9252e-24,8.8601e-24,8.7954e-24,8.7313e-24,8.6676e-24,8.6045e-24,8.5418e-24,8.4795e-24,8.4176e-24},
      std::array<float, 1001>{ 3.6588e-21,3.6479e-21,3.6368e-21,3.6254e-21,3.6137e-21,3.6019e-21,3.5898e-21,3.5774e-21,3.5649e-21,3.5522e-21,3.5392e-21,3.5261e-21,3.5128e-21,3.4993e-21,3.4857e-21,3.4719e-21,3.4579e-21,3.4438e-21,3.4297e-21,3.4154e-21,3.401e-21,3.3865e-21,3.3719e-21,3.3573e-21,3.3427e-21,3.3279e-21,3.3132e-21,3.2984e-21,3.2837e-21,3.2689e-21,3.2542e-21,3.2395e-21,3.2248e-21,3.2101e-21,3.1955e-21,3.181e-21,3.1665e-21,3.1522e-21,3.1379e-21,3.1237e-21,3.1095e-21,3.0955e-21,3.0817e-21,3.0679e-21,3.0543e-21,3.0407e-21,3.0274e-21,3.0141e-21,3.001e-21,2.9881e-21,2.9753e-21,2.9626e-21,2.9501e-21,2.9378e-21,2.9256e-21,2.9136e-21,2.9017e-21,2.89e-21,2.8785e-21,2.8671e-21,2.8559e-21,2.8448e-21,2.8339e-21,2.8232e-21,2.8127e-21,2.8023e-21,2.792e-21,2.782e-21,2.772e-21,2.7623e-21,2.7527e-21,2.7432e-21,2.7339e-21,2.7248e-21,2.7158e-21,2.7069e-21,2.6982e-21,2.6896e-21,2.6811e-21,2.6728e-21,2.6646e-21,2.6565e-21,2.6486e-21,2.6408e-21,2.6331e-21,2.6255e-21,2.618e-21,2.6106e-21,2.6034e-21,2.5962e-21,2.5891e-21,2.5822e-21,2.5753e-21,2.5685e-21,2.5617e-21,2.5551e-21,2.5485e-21,2.542e-21,2.5356e-21,2.5292e-21,2.5229e-21,2.5166e-21,2.5104e-21,2.5043e-21,2.4982e-21,2.4921e-21,2.4861e-21,2.4801e-21,2.4741e-21,2.4681e-21,2.4622e-21,2.4563e-21,2.4504e-21,2.4445e-21,2.4387e-21,2.4328e-21,2.4269e-21,2.4211e-21,2.4152e-21,2.4093e-21,2.4034e-21,2.3975e-21,2.3916e-21,2.3856e-21,2.3797e-21,2.3737e-21,2.3676e-21,2.3616e-21,2.3555e-21,2.3493e-21,2.3432e-21,2.337e-21,2.3307e-21,2.3244e-21,2.318e-21,2.3116e-21,2.3052e-21,2.2987e-21,2.2921e-21,2.2855e-21,2.2789e-21,2.2721e-21,2.2653e-21,2.2585e-21,2.2516e-21,2.2446e-21,2.2376e-21,2.2305e-21,2.2234e-21,2.2162e-21,2.2089e-21,2.2015e-21,2.1941e-21,2.1867e-21,2.1792e-21,2.1716e-21,2.1639e-21,2.1562e-21,2.1485e-21,2.1406e-21,2.1328e-21,2.1248e-21,2.1168e-21,2.1088e-21,2.1007e-21,2.0926e-21,2.0844e-21,2.0762e-21,2.0679e-21,2.0595e-21,2.0512e-21,2.0428e-21,2.0343e-21,2.0258e-21,2.0173e-21,2.0088e-21,2.0002e-21,1.9921e-21,1.9841e-21,1.976e-21,1.9679e-21,1.9598e-21,1.9516e-21,1.9435e-21,1.9353e-21,1.9271e-21,1.9189e-21,1.9107e-21,1.9024e-21,1.8942e-21,1.886e-21,1.8777e-21,1.8695e-21,1.8612e-21,1.853e-21,1.8448e-21,1.8365e-21,1.8283e-21,1.8201e-21,1.8119e-21,1.8037e-21,1.7955e-21,1.7873e-21,1.7792e-21,1.771e-21,1.7629e-21,1.7548e-21,1.7468e-21,1.7387e-21,1.7307e-21,1.7227e-21,1.7147e-21,1.7067e-21,1.6988e-21,1.6909e-21,1.6831e-21,1.6752e-21,1.6674e-21,1.6596e-21,1.6519e-21,1.6442e-21,1.6365e-21,1.6288e-21,1.6212e-21,1.6137e-21,1.6061e-21,1.5986e-21,1.5911e-21,1.5837e-21,1.5763e-21,1.5689e-21,1.5616e-21,1.5543e-21,1.547e-21,1.5398e-21,1.5326e-21,1.5255e-21,1.5184e-21,1.5113e-21,1.5042e-21,1.4972e-21,1.4903e-21,1.4833e-21,1.4765e-21,1.4696e-21,1.4628e-21,1.456e-21,1.4492e-21,1.4425e-21,1.4358e-21,1.4292e-21,1.4226e-21,1.416e-21,1.4094e-21,1.4029e-21,1.3965e-21,1.39e-21,1.3836e-21,1.3772e-21,1.3709e-21,1.3645e-21,1.3583e-21,1.352e-21,1.3458e-21,1.3396e-21,1.3334e-21,1.3273e-21,1.3212e-21,1.3151e-21,1.3091e-21,1.303e-21,1.297e-21,1.2911e-21,1.2851e-21,1.2792e-21,1.2733e-21,1.2675e-21,1.2616e-21,1.2558e-21,1.2501e-21,1.2443e-21,1.2386e-21,1.2329e-21,1.2272e-21,1.2215e-21,1.2159e-21,1.2103e-21,1.2047e-21,1.1991e-21,1.1935e-21,1.188e-21,1.1825e-21,1.177e-21,1.1715e-21,1.1661e-21,1.1607e-21,1.1553e-21,1.1499e-21,1.1445e-21,1.1392e-21,1.1339e-21,1.1285e-21,1.1237e-21,1.119e-21,1.1142e-21,1.1095e-21,1.1047e-21,1.1e-21,1.0953e-21,1.0906e-21,1.086e-21,1.0813e-21,1.0767e-21,1.0721e-21,1.0675e-21,1.0629e-21,1.0583e-21,1.0537e-21,1.0492e-21,1.0446e-21,1.0401e-21,1.0356e-21,1.0311e-21,1.0266e-21,1.0221e-21,1.0177e-21,1.0132e-21,1.0088e-21,1.0044e-21,9.9997e-22,9.9557e-22,9.9119e-22,9.8682e-22,9.8246e-22,9.7812e-22,9.7379e-22,9.6947e-22,9.6516e-22,9.6086e-22,9.5658e-22,9.523e-22,9.4804e-22,9.438e-22,9.3956e-22,9.3533e-22,9.3112e-22,9.2692e-22,9.2273e-22,9.1855e-22,9.1438e-22,9.1023e-22,9.0608e-22,9.0195e-22,8.9783e-22,8.9372e-22,8.8962e-22,8.8553e-22,8.8146e-22,8.7739e-22,8.7334e-22,8.693e-22,8.6526e-22,8.6124e-22,8.5723e-22,8.5324e-22,8.4925e-22,8.4527e-22,8.4131e-22,8.3735e-22,8.3341e-22,8.2948e-22,8.2555e-22,8.2164e-22,8.1774e-22,8.1385e-22,8.0997e-22,8.061e-22,8.0225e-22,7.984e-22,7.9456e-22,7.9074e-22,7.8692e-22,7.8312e-22,7.7933e-22,7.7554e-22,7.7177e-22,7.6801e-22,7.6426e-22,7.6052e-22,7.5679e-22,7.5306e-22,7.4936e-22,7.4566e-22,7.4197e-22,7.3829e-22,7.3462e-22,7.3096e-22,7.2732e-22,7.2371e-22,7.2053e-22,7.1736e-22,7.142e-22,7.1105e-22,7.0791e-22,7.0477e-22,7.0164e-22,6.9853e-22,6.9542e-22,6.9231e-22,6.8922e-22,6.8613e-22,6.8306e-22,6.7999e-22,6.7693e-22,6.7388e-22,6.7083e-22,6.678e-22,6.6477e-22,6.6175e-22,6.5874e-22,6.5574e-22,6.5274e-22,6.4976e-22,6.4678e-22,6.4381e-22,6.4085e-22,6.3789e-22,6.3495e-22,6.3201e-22,6.2908e-22,6.2616e-22,6.2325e-22,6.2034e-22,6.1745e-22,6.1456e-22,6.1168e-22,6.0881e-22,6.0594e-22,6.0309e-22,6.0024e-22,5.974e-22,5.9457e-22,5.9174e-22,5.8893e-22,5.8612e-22,5.8332e-22,5.8053e-22,5.7774e-22,5.7497e-22,5.722e-22,5.6944e-22,5.6669e-22,5.6394e-22,5.6121e-22,5.5848e-22,5.5576e-22,5.5304e-22,5.5034e-22,5.4764e-22,5.4495e-22,5.4227e-22,5.396e-22,5.3693e-22,5.3427e-22,5.3162e-22,5.2898e-22,5.2635e-22,5.2372e-22,5.211e-22,5.1849e-22,5.1589e-22,5.1329e-22,5.1071e-22,5.0813e-22,5.0555e-22,5.0299e-22,5.0043e-22,4.9788e-22,4.9571e-22,4.9359e-22,4.9148e-22,4.8938e-22,4.8728e-22,4.8518e-22,4.831e-22,4.8101e-22,4.7894e-22,4.7687e-22,4.748e-22,4.7275e-22,4.7069e-22,4.6865e-22,4.6661e-22,4.6457e-22,4.6254e-22,4.6052e-22,4.585e-22,4.5649e-22,4.5449e-22,4.5249e-22,4.5049e-22,4.485e-22,4.4652e-22,4.4455e-22,4.4258e-22,4.4061e-22,4.3865e-22,4.367e-22,4.3475e-22,4.3281e-22,4.3088e-22,4.2895e-22,4.2702e-22,4.2511e-22,4.232e-22,4.2129e-22,4.1939e-22,4.1749e-22,4.156e-22,4.1372e-22,4.1184e-22,4.0997e-22,4.0811e-22,4.0625e-22,4.0439e-22,4.0254e-22,4.007e-22,3.9886e-22,3.9703e-22,3.952e-22,3.9338e-22,3.9157e-22,3.8976e-22,3.8796e-22,3.8616e-22,3.8437e-22,3.8258e-22,3.808e-22,3.7902e-22,3.7725e-22,3.7549e-22,3.7373e-22,3.7198e-22,3.7023e-22,3.6849e-22,3.6712e-22,3.6578e-22,3.6445e-22,3.6312e-22,3.6179e-22,3.6047e-22,3.5915e-22,3.5783e-22,3.5652e-22,3.5521e-22,3.5391e-22,3.5261e-22,3.5132e-22,3.5002e-22,3.4874e-22,3.4745e-22,3.4617e-22,3.449e-22,3.4363e-22,3.4236e-22,3.4109e-22,3.3983e-22,3.3858e-22,3.3732e-22,3.3608e-22,3.3483e-22,3.3359e-22,3.3235e-22,3.3112e-22,3.2989e-22,3.2867e-22,3.2744e-22,3.2623e-22,3.2501e-22,3.238e-22,3.226e-22,3.2139e-22,3.2019e-22,3.19e-22,3.1781e-22,3.1662e-22,3.1544e-22,3.1426e-22,3.1308e-22,3.1191e-22,3.1074e-22,3.0958e-22,3.0841e-22,3.0726e-22,3.061e-22,3.0495e-22,3.0381e-22,3.0267e-22,3.0153e-22,3.0039e-22,2.9926e-22,2.9813e-22,2.9701e-22,2.9624e-22,2.955e-22,2.9476e-22,2.9402e-22,2.9329e-22,2.9255e-22,2.9182e-22,2.9109e-22,2.9037e-22,2.8964e-22,2.8892e-22,2.882e-22,2.8748e-22,2.8677e-22,2.8605e-22,2.8534e-22,2.8463e-22,2.8393e-22,2.8322e-22,2.8252e-22,2.8182e-22,2.8112e-22,2.8043e-22,2.7973e-22,2.7904e-22,2.7835e-22,2.7766e-22,2.7698e-22,2.763e-22,2.7562e-22,2.7494e-22,2.7426e-22,2.7359e-22,2.7292e-22,2.7225e-22,2.7158e-22,2.7091e-22,2.7025e-22,2.6959e-22,2.6893e-22,2.6827e-22,2.6762e-22,2.6696e-22,2.6631e-22,2.6566e-22,2.6502e-22,2.6437e-22,2.6373e-22,2.6309e-22,2.6245e-22,2.6182e-22,2.6219e-22,2.6283e-22,2.6346e-22,2.6409e-22,2.6472e-22,2.6534e-22,2.6597e-22,2.6659e-22,2.6721e-22,2.6782e-22,2.6843e-22,2.6905e-22,2.6965e-22,2.7026e-22,2.7087e-22,2.7147e-22,2.7207e-22,2.7266e-22,2.7326e-22,2.7385e-22,2.7444e-22,2.7503e-22,2.7562e-22,2.762e-22,2.7679e-22,2.7737e-22,2.7795e-22,2.7852e-22,2.791e-22,2.7967e-22,2.8024e-22,2.8081e-22,2.8137e-22,2.8194e-22,2.825e-22,2.8306e-22,2.8362e-22,2.8417e-22,2.8473e-22,2.8528e-22,2.8583e-22,2.8638e-22,2.8693e-22,2.8747e-22,2.8801e-22,2.8849e-22,2.874e-22,2.863e-22,2.8522e-22,2.8413e-22,2.8305e-22,2.8198e-22,2.809e-22,2.7983e-22,2.7877e-22,2.7771e-22,2.7665e-22,2.7559e-22,2.7454e-22,2.7349e-22,2.7245e-22,2.7141e-22,2.7037e-22,2.6933e-22,2.683e-22,2.6728e-22,2.6625e-22,2.6523e-22,2.6422e-22,2.632e-22,2.6219e-22,2.6119e-22,2.6018e-22,2.5918e-22,2.5819e-22,2.5719e-22,2.562e-22,2.5522e-22,2.5423e-22,2.5325e-22,2.5227e-22,2.513e-22,2.5033e-22,2.4936e-22,2.484e-22,2.4744e-22,2.4648e-22,2.4634e-22,2.4664e-22,2.4694e-22,2.4724e-22,2.4754e-22,2.4784e-22,2.4814e-22,2.4844e-22,2.4873e-22,2.4903e-22,2.4933e-22,2.4962e-22,2.4991e-22,2.502e-22,2.505e-22,2.5079e-22,2.5108e-22,2.5137e-22,2.5165e-22,2.5194e-22,2.5223e-22,2.5251e-22,2.528e-22,2.5308e-22,2.5337e-22,2.5365e-22,2.5393e-22,2.5421e-22,2.5449e-22,2.5477e-22,2.5505e-22,2.5533e-22,2.5561e-22,2.5589e-22,2.5616e-22,2.5644e-22,2.5671e-22,2.5699e-22,2.5715e-22,2.573e-22,2.5745e-22,2.5759e-22,2.5774e-22,2.5789e-22,2.5803e-22,2.5818e-22,2.5832e-22,2.5847e-22,2.5861e-22,2.5876e-22,2.589e-22,2.5904e-22,2.5919e-22,2.5933e-22,2.5947e-22,2.5962e-22,2.5976e-22,2.599e-22,2.6004e-22,2.6018e-22,2.6032e-22,2.6046e-22,2.6061e-22,2.6075e-22,2.6089e-22,2.6103e-22,2.6117e-22,2.613e-22,2.6144e-22,2.6158e-22,2.6172e-22,2.6186e-22,2.6198e-22,2.619e-22,2.6183e-22,2.6175e-22,2.6168e-22,2.6161e-22,2.6154e-22,2.6147e-22,2.6139e-22,2.6132e-22,2.6125e-22,2.6118e-22,2.6111e-22,2.6104e-22,2.6097e-22,2.609e-22,2.6083e-22,2.6077e-22,2.607e-22,2.6063e-22,2.6056e-22,2.605e-22,2.6043e-22,2.6036e-22,2.603e-22,2.6023e-22,2.6017e-22,2.601e-22,2.6004e-22,2.5998e-22,2.5991e-22,2.5985e-22,2.5979e-22,2.602e-22,2.6066e-22,2.6113e-22,2.6159e-22,2.6204e-22,2.625e-22,2.6296e-22,2.6341e-22,2.6387e-22,2.6432e-22,2.6477e-22,2.6522e-22,2.6567e-22,2.6612e-22,2.6656e-22,2.6701e-22,2.6746e-22,2.679e-22,2.6834e-22,2.6878e-22,2.6922e-22,2.6966e-22,2.701e-22,2.7054e-22,2.7097e-22,2.7141e-22,2.7184e-22,2.7227e-22,2.7271e-22,2.7314e-22,2.7818e-22,2.8345e-22,2.887e-22,2.9395e-22,2.9918e-22,3.044e-22,3.0961e-22,3.148e-22,3.1999e-22,3.2516e-22,3.3032e-22,3.3547e-22,3.4061e-22,3.4573e-22,3.5085e-22,3.5595e-22,3.6104e-22,3.6612e-22,3.7119e-22,3.7625e-22,3.813e-22,3.8633e-22,3.9136e-22,3.9637e-22,4.0137e-22,4.0636e-22,4.1134e-22,4.1631e-22,4.1536e-22,4.1391e-22,4.1247e-22,4.1103e-22,4.096e-22,4.0817e-22,4.0674e-22,4.0532e-22,4.039e-22,4.0249e-22,4.0108e-22,3.9968e-22,3.9827e-22,3.9688e-22,3.9548e-22,3.9409e-22,3.9271e-22,3.9133e-22,3.8995e-22,3.8857e-22,3.872e-22,3.8584e-22,3.8448e-22,3.8312e-22,3.8176e-22,3.8041e-22,3.7906e-22,3.7772e-22,3.7638e-22,3.7505e-22,3.7371e-22,3.7239e-22,3.7106e-22,3.6974e-22,3.6842e-22,3.6711e-22,3.658e-22,3.645e-22,3.6319e-22,3.6189e-22,3.606e-22,3.5931e-22,3.5802e-22,3.5674e-22,3.5546e-22,3.5418e-22,3.5291e-22,3.5164e-22,3.5037e-22,3.4911e-22,3.4785e-22,3.4868e-22,3.5014e-22,3.516e-22,3.5305e-22,3.5451e-22,3.5596e-22,3.5741e-22,3.5886e-22,3.603e-22,3.6175e-22,3.6319e-22,3.6464e-22,3.6608e-22,3.6752e-22,3.6896e-22,3.704e-22,3.7183e-22,3.7327e-22,3.7471e-22,3.7614e-22,3.7757e-22,3.79e-22,3.8043e-22,3.8186e-22,3.8329e-22,3.8472e-22,3.8614e-22,3.8757e-22,3.8899e-22,3.9042e-22,3.9184e-22,3.9326e-22,3.9468e-22,3.961e-22,3.9752e-22,3.9894e-22,4.0036e-22,4.0177e-22,4.0319e-22,4.046e-22,4.0602e-22,4.0743e-22,4.0884e-22,4.1026e-22,4.1167e-22,4.131e-22}};

    const static float base_density = 8.544e+18f; //km^-3

    const static float Hp = 730e-3f; // Continental value

    const static float pb = 2e6f; // background value in km^-3

    template <typename Float, typename UInt32, typename Mask, typename Spectrum, typename Wavelength>
    class RuralAerosol : public GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength> {
        const float pb_bd;

    public:
        RuralAerosol() : GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength>(tabulatedValues),
                         pb_bd(pb / base_density) {}

        Float get_density(const Float &z) const override {
            return Float(base_density) * (enoki::exp(-z / Float(Hp)) + Float(pb_bd));
        }

        [[nodiscard]] float get_density_float(const float &z) const override {
            return base_density * (enoki::exp(-z / Hp) + pb_bd);
        }
    };
}
#endif //_AEROSOL_RURAL