#ifndef _AEROSOL_URBAN
#define _AEROSOL_URBAN

#include <array>
#include <mitsuba/atmosphere/GlobalAtmosphericAerosol.h>

namespace urbanAerosol {

	// row[0] = wavelength
	// row[1] = cross section absorption coefficient
	// row[2] = cross section scattering coefficient
    const static std::array<std::array<float, 1001>, 3> tabulatedValues =
	{ std::array<float, 1001>{ 100,100.23,100.46,100.69,100.93,101.16,101.39,101.62,101.86,102.09,102.33,102.57,102.8,103.04,103.28,103.51,103.75,103.99,104.23,104.47,104.71,104.95,105.2,105.44,105.68,105.93,106.17,106.41,106.66,106.91,107.15,107.4,107.65,107.89,108.14,108.39,108.64,108.89,109.14,109.4,109.65,109.9,110.15,110.41,110.66,110.92,111.17,111.43,111.69,111.94,112.2,112.46,112.72,112.98,113.24,113.5,113.76,114.03,114.29,114.55,114.82,115.08,115.35,115.61,115.88,116.14,116.41,116.68,116.95,117.22,117.49,117.76,118.03,118.3,118.58,118.85,119.12,119.4,119.67,119.95,120.23,120.5,120.78,121.06,121.34,121.62,121.9,122.18,122.46,122.74,123.03,123.31,123.59,123.88,124.17,124.45,124.74,125.03,125.31,125.6,125.89,126.18,126.47,126.77,127.06,127.35,127.64,127.94,128.23,128.53,128.82,129.12,129.42,129.72,130.02,130.32,130.62,130.92,131.22,131.52,131.83,132.13,132.43,132.74,133.05,133.35,133.66,133.97,134.28,134.59,134.9,135.21,135.52,135.83,136.14,136.46,136.77,137.09,137.4,137.72,138.04,138.36,138.68,139,139.32,139.64,139.96,140.28,140.6,140.93,141.25,141.58,141.91,142.23,142.56,142.89,143.22,143.55,143.88,144.21,144.54,144.88,145.21,145.55,145.88,146.22,146.55,146.89,147.23,147.57,147.91,148.25,148.59,148.94,149.28,149.62,149.97,150.31,150.66,151.01,151.36,151.71,152.05,152.41,152.76,153.11,153.46,153.82,154.17,154.53,154.88,155.24,155.6,155.96,156.31,156.68,157.04,157.4,157.76,158.12,158.49,158.85,159.22,159.59,159.96,160.32,160.69,161.06,161.44,161.81,162.18,162.55,162.93,163.31,163.68,164.06,164.44,164.82,165.2,165.58,165.96,166.34,166.72,167.11,167.49,167.88,168.27,168.66,169.04,169.43,169.82,170.22,170.61,171,171.4,171.79,172.19,172.58,172.98,173.38,173.78,174.18,174.58,174.98,175.39,175.79,176.2,176.6,177.01,177.42,177.83,178.24,178.65,179.06,179.47,179.89,180.3,180.72,181.13,181.55,181.97,182.39,182.81,183.23,183.65,184.08,184.5,184.93,185.35,185.78,186.21,186.64,187.07,187.5,187.93,188.36,188.8,189.23,189.67,190.11,190.55,190.99,191.43,191.87,192.31,192.75,193.2,193.64,194.09,194.54,194.98,195.43,195.88,196.34,196.79,197.24,197.7,198.15,198.61,199.07,199.53,199.99,200.45,200.91,201.37,201.84,202.3,202.77,203.24,203.7,204.17,204.64,205.12,205.59,206.06,206.54,207.01,207.49,207.97,208.45,208.93,209.41,209.89,210.38,210.86,211.35,211.84,212.32,212.81,213.3,213.8,214.29,214.78,215.28,215.77,216.27,216.77,217.27,217.77,218.27,218.78,219.28,219.79,220.29,220.8,221.31,221.82,222.33,222.84,223.36,223.87,224.39,224.91,225.42,225.94,226.46,226.99,227.51,228.03,228.56,229.09,229.61,230.14,230.67,231.21,231.74,232.27,232.81,233.35,233.88,234.42,234.96,235.5,236.05,236.59,237.14,237.68,238.23,238.78,239.33,239.88,240.44,240.99,241.55,242.1,242.66,243.22,243.78,244.34,244.91,245.47,246.04,246.6,247.17,247.74,248.31,248.89,249.46,250.03,250.61,251.19,251.77,252.35,252.93,253.51,254.1,254.68,255.27,255.86,256.45,257.04,257.63,258.23,258.82,259.42,260.02,260.62,261.22,261.82,262.42,263.03,263.63,264.24,264.85,265.46,266.07,266.69,267.3,267.92,268.53,269.15,269.77,270.4,271.02,271.64,272.27,272.9,273.53,274.16,274.79,275.42,276.06,276.69,277.33,277.97,278.61,279.25,279.9,280.54,281.19,281.84,282.49,283.14,283.79,284.45,285.1,285.76,286.42,287.08,287.74,288.4,289.07,289.73,290.4,291.07,291.74,292.42,293.09,293.76,294.44,295.12,295.8,296.48,297.17,297.85,298.54,299.23,299.92,300.61,301.3,302,302.69,303.39,304.09,304.79,305.49,306.2,306.9,307.61,308.32,309.03,309.74,310.46,311.17,311.89,312.61,313.33,314.05,314.77,315.5,316.23,316.96,317.69,318.42,319.15,319.89,320.63,321.37,322.11,322.85,323.59,324.34,325.09,325.84,326.59,327.34,328.1,328.85,329.61,330.37,331.13,331.89,332.66,333.43,334.19,334.97,335.74,336.51,337.29,338.06,338.84,339.63,340.41,341.19,341.98,342.77,343.56,344.35,345.14,345.94,346.74,347.54,348.34,349.14,349.95,350.75,351.56,352.37,353.18,354,354.81,355.63,356.45,357.27,358.1,358.92,359.75,360.58,361.41,362.24,363.08,363.92,364.75,365.59,366.44,367.28,368.13,368.98,369.83,370.68,371.54,372.39,373.25,374.11,374.97,375.84,376.7,377.57,378.44,379.32,380.19,381.07,381.94,382.82,383.71,384.59,385.48,386.37,387.26,388.15,389.05,389.94,390.84,391.74,392.64,393.55,394.46,395.37,396.28,397.19,398.11,399.02,399.94,400.87,401.79,402.72,403.65,404.58,405.51,406.44,407.38,408.32,409.26,410.2,411.15,412.1,413.05,414,414.95,415.91,416.87,417.83,418.79,419.76,420.73,421.7,422.67,423.64,424.62,425.6,426.58,427.56,428.55,429.54,430.53,431.52,432.51,433.51,434.51,435.51,436.52,437.52,438.53,439.54,440.55,441.57,442.59,443.61,444.63,445.66,446.68,447.71,448.75,449.78,450.82,451.86,452.9,453.94,454.99,456.04,457.09,458.14,459.2,460.26,461.32,462.38,463.45,464.52,465.59,466.66,467.74,468.81,469.89,470.98,472.06,473.15,474.24,475.34,476.43,477.53,478.63,479.73,480.84,481.95,483.06,484.17,485.29,486.41,487.53,488.65,489.78,490.91,492.04,493.17,494.31,495.45,496.59,497.74,498.88,500.03,501.19,502.34,503.5,504.66,505.82,506.99,508.16,509.33,510.51,511.68,512.86,514.04,515.23,516.42,517.61,518.8,520,521.19,522.4,523.6,524.81,526.02,527.23,528.45,529.66,530.88,532.11,533.33,534.56,535.8,537.03,538.27,539.51,540.75,542,543.25,544.5,545.76,547.02,548.28,549.54,550.81,552.08,553.35,554.63,555.9,557.19,558.47,559.76,561.05,562.34,563.64,564.94,566.24,567.54,568.85,570.16,571.48,572.8,574.12,575.44,576.77,578.1,579.43,580.76,582.1,583.45,584.79,586.14,587.49,588.84,590.2,591.56,592.93,594.29,595.66,597.04,598.41,599.79,601.17,602.56,603.95,605.34,606.74,608.13,609.54,610.94,612.35,613.76,615.18,616.6,618.02,619.44,620.87,622.3,623.73,625.17,626.61,628.06,629.51,630.96,632.41,633.87,635.33,636.8,638.26,639.73,641.21,642.69,644.17,645.65,647.14,648.63,650.13,651.63,653.13,654.64,656.15,657.66,659.17,660.69,662.22,663.74,665.27,666.81,668.34,669.88,671.43,672.98,674.53,676.08,677.64,679.2,680.77,682.34,683.91,685.49,687.07,688.65,690.24,691.83,693.43,695.02,696.63,698.23,699.84,701.46,703.07,704.69,706.32,707.95,709.58,711.21,712.85,714.5,716.14,717.79,719.45,721.11,722.77,724.44,726.11,727.78,729.46,731.14,732.82,734.51,736.21,737.9,739.61,741.31,743.02,744.73,746.45,748.17,749.89,751.62,753.36,755.09,756.83,758.58,760.33,762.08,763.84,765.6,767.36,769.13,770.9,772.68,774.46,776.25,778.04,779.83,781.63,783.43,785.24,787.05,788.86,790.68,792.5,794.33,796.16,797.99,799.83,801.68,803.53,805.38,807.24,809.1,810.96,812.83,814.7,816.58,818.46,820.35,822.24,824.14,826.04,827.94,829.85,831.76,833.68,835.6,837.53,839.46,841.4,843.33,845.28,847.23,849.18,851.14,853.1,855.07,857.04,859.01,860.99,862.98,864.97,866.96,868.96,870.96,872.97,874.98,877,879.02,881.05,883.08,885.12,887.16,889.2,891.25,893.31,895.36,897.43,899.5,901.57,903.65,905.73,907.82,909.91,912.01,914.11,916.22,918.33,920.45,922.57,924.7,926.83,928.97,931.11,933.25,935.41,937.56,939.72,941.89,944.06,946.24,948.42,950.6,952.8,954.99,957.19,959.4,961.61,963.83,966.05,968.28,970.51,972.75,974.99,977.24,979.49,981.75,984.01,986.28,988.55,990.83,993.12,995.41,997.7,1000},
      std::array<float, 1001>{ 1.8996e-21,1.8909e-21,1.8822e-21,1.8735e-21,1.8649e-21,1.8563e-21,1.8478e-21,1.8394e-21,1.8309e-21,1.8226e-21,1.8142e-21,1.8059e-21,1.7977e-21,1.7895e-21,1.7813e-21,1.7732e-21,1.7651e-21,1.7571e-21,1.7491e-21,1.7411e-21,1.7332e-21,1.7253e-21,1.7174e-21,1.7096e-21,1.7018e-21,1.6941e-21,1.6864e-21,1.6787e-21,1.671e-21,1.6634e-21,1.6558e-21,1.6483e-21,1.6407e-21,1.6332e-21,1.6258e-21,1.6183e-21,1.6109e-21,1.6035e-21,1.5961e-21,1.5888e-21,1.5815e-21,1.5742e-21,1.5669e-21,1.5597e-21,1.5525e-21,1.5453e-21,1.5381e-21,1.531e-21,1.5239e-21,1.5168e-21,1.5097e-21,1.5026e-21,1.4956e-21,1.4886e-21,1.4816e-21,1.4746e-21,1.4676e-21,1.4607e-21,1.4537e-21,1.4468e-21,1.4399e-21,1.4331e-21,1.4262e-21,1.4194e-21,1.4126e-21,1.4057e-21,1.399e-21,1.3922e-21,1.3854e-21,1.3787e-21,1.372e-21,1.3653e-21,1.3586e-21,1.3519e-21,1.3452e-21,1.3386e-21,1.3319e-21,1.3253e-21,1.3187e-21,1.3121e-21,1.3055e-21,1.299e-21,1.2924e-21,1.2859e-21,1.2794e-21,1.2729e-21,1.2664e-21,1.2599e-21,1.2534e-21,1.247e-21,1.2405e-21,1.2341e-21,1.2277e-21,1.2213e-21,1.2149e-21,1.2085e-21,1.2022e-21,1.1959e-21,1.1895e-21,1.1832e-21,1.1769e-21,1.1706e-21,1.1643e-21,1.1581e-21,1.1518e-21,1.1456e-21,1.1394e-21,1.1332e-21,1.127e-21,1.1208e-21,1.1146e-21,1.1085e-21,1.1024e-21,1.0962e-21,1.0901e-21,1.084e-21,1.078e-21,1.0719e-21,1.0658e-21,1.0598e-21,1.0538e-21,1.0478e-21,1.0418e-21,1.0358e-21,1.0299e-21,1.0239e-21,1.018e-21,1.0121e-21,1.0062e-21,1.0003e-21,9.9443e-22,9.8857e-22,9.8274e-22,9.7692e-22,9.7111e-22,9.6533e-22,9.5956e-22,9.538e-22,9.4806e-22,9.4234e-22,9.3663e-22,9.3094e-22,9.2527e-22,9.1962e-22,9.1398e-22,9.0836e-22,9.0275e-22,8.9717e-22,8.916e-22,8.8605e-22,8.8051e-22,8.75e-22,8.695e-22,8.6402e-22,8.5856e-22,8.5311e-22,8.4769e-22,8.4228e-22,8.3689e-22,8.3152e-22,8.2617e-22,8.2084e-22,8.1553e-22,8.1023e-22,8.0496e-22,7.997e-22,7.9447e-22,7.8925e-22,7.8405e-22,7.7887e-22,7.7371e-22,7.6857e-22,7.6346e-22,7.5836e-22,7.5328e-22,7.4822e-22,7.4318e-22,7.3816e-22,7.3316e-22,7.2818e-22,7.2322e-22,7.1828e-22,7.1336e-22,7.0847e-22,7.0359e-22,6.9873e-22,6.939e-22,6.8909e-22,6.8429e-22,6.7952e-22,6.7477e-22,6.7004e-22,6.6534e-22,6.6065e-22,6.5599e-22,6.5134e-22,6.4672e-22,6.4212e-22,6.3755e-22,6.3299e-22,6.2846e-22,6.2394e-22,6.1945e-22,6.1498e-22,6.1054e-22,6.0611e-22,6.0171e-22,5.9733e-22,5.9297e-22,5.8863e-22,5.8432e-22,5.8003e-22,5.7576e-22,5.7151e-22,5.6728e-22,5.6308e-22,5.589e-22,5.5474e-22,5.506e-22,5.4649e-22,5.4239e-22,5.3832e-22,5.3428e-22,5.3025e-22,5.2625e-22,5.2226e-22,5.1831e-22,5.1437e-22,5.1045e-22,5.0656e-22,5.0269e-22,4.9884e-22,4.9502e-22,4.9121e-22,4.8743e-22,4.8367e-22,4.7994e-22,4.7622e-22,4.7253e-22,4.6886e-22,4.6521e-22,4.6158e-22,4.5798e-22,4.5439e-22,4.5083e-22,4.4729e-22,4.4377e-22,4.4028e-22,4.368e-22,4.3335e-22,4.2992e-22,4.2651e-22,4.2312e-22,4.1976e-22,4.1641e-22,4.1309e-22,4.0979e-22,4.0651e-22,4.0325e-22,4.0001e-22,3.9679e-22,3.936e-22,3.9042e-22,3.8727e-22,3.8413e-22,3.8102e-22,3.7793e-22,3.7486e-22,3.7181e-22,3.6878e-22,3.6577e-22,3.6278e-22,3.5981e-22,3.5687e-22,3.5394e-22,3.5103e-22,3.4814e-22,3.4528e-22,3.4243e-22,3.396e-22,3.3679e-22,3.3401e-22,3.3124e-22,3.2849e-22,3.2576e-22,3.2305e-22,3.2036e-22,3.1769e-22,3.1503e-22,3.124e-22,3.0979e-22,3.0719e-22,3.0462e-22,3.0206e-22,2.9952e-22,2.97e-22,2.945e-22,2.9201e-22,2.8955e-22,2.871e-22,2.8467e-22,2.8226e-22,2.7987e-22,2.7749e-22,2.7513e-22,2.7279e-22,2.7047e-22,2.6817e-22,2.6588e-22,2.6361e-22,2.6136e-22,2.5912e-22,2.569e-22,2.547e-22,2.5251e-22,2.5035e-22,2.4819e-22,2.4606e-22,2.4394e-22,2.4184e-22,2.3975e-22,2.3768e-22,2.3563e-22,2.3359e-22,2.3157e-22,2.2957e-22,2.2758e-22,2.256e-22,2.2364e-22,2.217e-22,2.1977e-22,2.1786e-22,2.1596e-22,2.1408e-22,2.1222e-22,2.1036e-22,2.0853e-22,2.067e-22,2.049e-22,2.031e-22,2.0132e-22,1.9956e-22,1.9781e-22,1.9607e-22,1.9435e-22,1.9264e-22,1.9095e-22,1.8927e-22,1.876e-22,1.8595e-22,1.8431e-22,1.8268e-22,1.8107e-22,1.7947e-22,1.7788e-22,1.7631e-22,1.7475e-22,1.732e-22,1.7167e-22,1.7015e-22,1.6864e-22,1.6714e-22,1.6566e-22,1.6418e-22,1.6272e-22,1.6128e-22,1.5984e-22,1.5842e-22,1.5701e-22,1.5561e-22,1.5422e-22,1.5284e-22,1.5148e-22,1.5013e-22,1.4878e-22,1.4745e-22,1.4613e-22,1.4483e-22,1.4353e-22,1.4224e-22,1.4097e-22,1.3971e-22,1.3845e-22,1.3721e-22,1.3598e-22,1.3476e-22,1.3355e-22,1.3235e-22,1.3116e-22,1.2998e-22,1.2881e-22,1.2765e-22,1.265e-22,1.2536e-22,1.2423e-22,1.2311e-22,1.2199e-22,1.2089e-22,1.198e-22,1.1872e-22,1.1765e-22,1.1658e-22,1.1553e-22,1.1448e-22,1.1345e-22,1.1242e-22,1.114e-22,1.1039e-22,1.0939e-22,1.084e-22,1.0742e-22,1.0644e-22,1.0548e-22,1.0452e-22,1.0357e-22,1.0263e-22,1.017e-22,1.0078e-22,9.9859e-23,9.8951e-23,9.8051e-23,9.7159e-23,9.6275e-23,9.5399e-23,9.453e-23,9.3669e-23,9.2815e-23,9.1969e-23,9.1131e-23,9.03e-23,8.9476e-23,8.866e-23,8.7851e-23,8.7049e-23,8.6254e-23,8.5467e-23,8.4686e-23,8.3912e-23,8.3145e-23,8.2385e-23,8.1632e-23,8.0885e-23,8.0146e-23,7.9412e-23,7.8686e-23,7.7965e-23,7.7251e-23,7.6544e-23,7.5843e-23,7.5148e-23,7.446e-23,7.3777e-23,7.3101e-23,7.2431e-23,7.1766e-23,7.1108e-23,7.0456e-23,6.9809e-23,6.9168e-23,6.8534e-23,6.7904e-23,6.7281e-23,6.6663e-23,6.605e-23,6.5444e-23,6.4842e-23,6.4246e-23,6.3656e-23,6.307e-23,6.249e-23,6.1916e-23,6.1346e-23,6.0782e-23,6.0223e-23,5.9668e-23,5.9119e-23,5.8575e-23,5.8036e-23,5.7501e-23,5.6972e-23,5.6447e-23,5.5927e-23,5.5412e-23,5.4901e-23,5.4395e-23,5.3894e-23,5.3397e-23,5.2904e-23,5.2417e-23,5.1933e-23,5.1454e-23,5.0979e-23,5.0509e-23,5.0043e-23,4.9581e-23,4.9124e-23,4.867e-23,4.8221e-23,4.7776e-23,4.7335e-23,4.6897e-23,4.6464e-23,4.6035e-23,4.561e-23,4.5189e-23,4.4771e-23,4.4357e-23,4.3947e-23,4.3541e-23,4.3139e-23,4.274e-23,4.2345e-23,4.1953e-23,4.1565e-23,4.1181e-23,4.08e-23,4.0423e-23,4.0049e-23,3.9678e-23,3.9311e-23,3.8947e-23,3.8587e-23,3.823e-23,3.7876e-23,3.7526e-23,3.7178e-23,3.6834e-23,3.6493e-23,3.6155e-23,3.582e-23,3.5489e-23,3.516e-23,3.4834e-23,3.4512e-23,3.4192e-23,3.3875e-23,3.3561e-23,3.325e-23,3.2942e-23,3.2637e-23,3.2335e-23,3.2035e-23,3.1738e-23,3.1444e-23,3.1153e-23,3.0864e-23,3.0578e-23,3.0294e-23,3.0013e-23,2.9735e-23,2.9459e-23,2.9186e-23,2.8915e-23,2.8647e-23,2.8381e-23,2.8118e-23,2.7857e-23,2.7599e-23,2.7343e-23,2.7089e-23,2.6838e-23,2.6589e-23,2.6342e-23,2.6098e-23,2.5856e-23,2.5616e-23,2.5378e-23,2.5142e-23,2.4909e-23,2.4678e-23,2.4449e-23,2.4222e-23,2.3997e-23,2.3774e-23,2.3554e-23,2.3335e-23,2.3118e-23,2.2904e-23,2.2691e-23,2.248e-23,2.2271e-23,2.2065e-23,2.186e-23,2.1657e-23,2.1456e-23,2.1256e-23,2.1059e-23,2.0863e-23,2.0669e-23,2.0477e-23,2.0287e-23,2.0099e-23,1.9912e-23,1.9727e-23,1.9544e-23,1.9362e-23,1.9182e-23,1.9004e-23,1.8827e-23,1.8653e-23,1.8479e-23,1.8308e-23,1.8137e-23,1.7969e-23,1.7802e-23,1.7636e-23,1.7472e-23,1.731e-23,1.7149e-23,1.699e-23,1.6832e-23,1.6675e-23,1.652e-23,1.6367e-23,1.6215e-23,1.6064e-23,1.5915e-23,1.5767e-23,1.562e-23,1.5475e-23,1.5331e-23,1.5188e-23,1.5047e-23,1.4907e-23,1.4769e-23,1.4631e-23,1.4495e-23,1.4361e-23,1.4227e-23,1.4095e-23,1.3964e-23,1.3834e-23,1.3705e-23,1.3578e-23,1.3452e-23,1.3326e-23,1.3202e-23,1.308e-23,1.2958e-23,1.2838e-23,1.2718e-23,1.26e-23,1.2483e-23,1.2367e-23,1.2252e-23,1.2138e-23,1.2025e-23,1.1913e-23,1.1802e-23,1.1692e-23,1.1583e-23,1.1476e-23,1.1369e-23,1.1263e-23,1.1158e-23,1.1055e-23,1.0952e-23,1.085e-23,1.0749e-23,1.0649e-23,1.055e-23,1.0452e-23,1.0355e-23,1.0258e-23,1.0163e-23,1.0068e-23,9.9747e-24,9.8819e-24,9.79e-24,9.6989e-24,9.6087e-24,9.5193e-24,9.4308e-24,9.343e-24,9.2561e-24,9.17e-24,9.0847e-24,9.0002e-24,8.9165e-24,8.8336e-24,8.7514e-24,8.67e-24,8.5893e-24,8.5094e-24,8.4303e-24,8.3519e-24,8.2742e-24,8.1972e-24,8.121e-24,8.0454e-24,7.9706e-24,7.8964e-24,7.823e-24,7.7502e-24,7.6781e-24,7.6067e-24,7.536e-24,7.4659e-24,7.3964e-24,7.3276e-24,7.2588e-24,7.1906e-24,7.123e-24,7.0561e-24,6.9899e-24,6.9242e-24,6.8591e-24,6.7947e-24,6.7309e-24,6.6676e-24,6.605e-24,6.5429e-24,6.4815e-24,6.4206e-24,6.3602e-24,6.3005e-24,6.2413e-24,6.1826e-24,6.1245e-24,6.067e-24,6.01e-24,5.9535e-24,5.8975e-24,5.8421e-24,5.7872e-24,5.7328e-24,5.6789e-24,5.6255e-24,5.5727e-24,5.5203e-24,5.4684e-24,5.417e-24,5.3661e-24,5.3156e-24,5.2657e-24,5.2162e-24,5.1671e-24,5.1185e-24,5.0704e-24,5.0227e-24,4.9755e-24,4.9291e-24,4.8832e-24,4.8378e-24,4.7927e-24,4.7481e-24,4.704e-24,4.6602e-24,4.6168e-24,4.5739e-24,4.5313e-24,4.4891e-24,4.4474e-24,4.406e-24,4.365e-24,4.3244e-24,4.2841e-24,4.2443e-24,4.2048e-24,4.1657e-24,4.1269e-24,4.0885e-24,4.0505e-24,4.0128e-24,3.9754e-24,3.9385e-24,3.9018e-24,3.8655e-24,3.8295e-24,3.7939e-24,3.7586e-24,3.7236e-24,3.689e-24,3.6547e-24,3.6207e-24,3.587e-24,3.5536e-24,3.5205e-24,3.4878e-24,3.4556e-24,3.4238e-24,3.3923e-24,3.3611e-24,3.3301e-24,3.2995e-24,3.2691e-24,3.239e-24,3.2092e-24,3.1797e-24,3.1504e-24,3.1214e-24,3.0927e-24,3.0642e-24,3.036e-24,3.0081e-24,2.9804e-24,2.953e-24,2.9258e-24,2.8989e-24,2.8722e-24,2.8458e-24,2.8196e-24,2.7937e-24,2.768e-24,2.7425e-24,2.7173e-24,2.6923e-24,2.6675e-24,2.643e-24,2.6186e-24,2.5946e-24,2.5707e-24,2.547e-24,2.5236e-24,2.5001e-24,2.4768e-24,2.4538e-24,2.4309e-24,2.4083e-24,2.3859e-24,2.3637e-24,2.3417e-24,2.3199e-24,2.2983e-24,2.2769e-24,2.2557e-24,2.2347e-24,2.2139e-24,2.1933e-24,2.1729e-24,2.1526e-24,2.1326e-24,2.1128e-24,2.0931e-24,2.0736e-24,2.0543e-24,2.0352e-24,2.0162e-24,1.9975e-24,1.9789e-24,1.9604e-24,1.9422e-24,1.9241e-24,1.9062e-24,1.8885e-24,1.8709e-24,1.8537e-24,1.8366e-24,1.8198e-24,1.8031e-24,1.7865e-24,1.7701e-24,1.7538e-24,1.7377e-24,1.7218e-24,1.7059e-24,1.6903e-24,1.6748e-24,1.6594e-24,1.6441e-24,1.629e-24,1.6141e-24,1.5993e-24,1.5846e-24,1.57e-24,1.5556e-24,1.5413e-24,1.5272e-24,1.5132e-24,1.4993e-24,1.4855e-24,1.4719e-24,1.4584e-24,1.445e-24,1.4317e-24,1.4186e-24,1.4051e-24,1.3918e-24,1.3786e-24,1.3655e-24,1.3525e-24,1.3397e-24,1.327e-24,1.3144e-24,1.302e-24,1.2896e-24,1.2774e-24,1.2653e-24,1.2534e-24,1.2415e-24,1.2298e-24,1.2182e-24,1.2067e-24,1.1953e-24,1.1841e-24,1.1729e-24,1.1619e-24,1.1509e-24,1.1401e-24,1.1294e-24,1.1188e-24,1.1083e-24,1.0978e-24,1.0875e-24,1.0776e-24,1.0677e-24,1.058e-24,1.0483e-24,1.0387e-24,1.0292e-24,1.0198e-24,1.0104e-24,1.0012e-24,9.9205e-25,9.8298e-25,9.7399e-25,9.6508e-25,9.5625e-25,9.4751e-25,9.3884e-25,9.3026e-25,9.2175e-25,9.1332e-25,9.0497e-25,8.967e-25,8.885e-25,8.8038e-25,8.7233e-25,8.6435e-25,8.5645e-25,8.4862e-25,8.4086e-25,8.3318e-25,8.2556e-25,8.1802e-25,8.1054e-25,8.0313e-25,7.9579e-25,7.8852e-25,7.8131e-25,7.7417e-25,7.6709e-25,7.6008e-25,7.5314e-25,7.4626e-25,7.3944e-25,7.3268e-25,7.2598e-25,7.1935e-25,7.1278e-25,7.0627e-25,6.9981e-25,6.9342e-25,6.8708e-25,6.8081e-25,6.7798e-25,6.7617e-25,6.7435e-25,6.7253e-25,6.707e-25,6.6886e-25,6.6702e-25,6.6517e-25,6.6332e-25,6.6145e-25,6.5958e-25,6.577e-25,6.5582e-25,6.5392e-25,6.5202e-25,6.5011e-25,6.4819e-25,6.4626e-25,6.4432e-25,6.4237e-25,6.4042e-25,6.3845e-25,6.3648e-25,6.3449e-25,6.325e-25,6.305e-25,6.2849e-25,6.2646e-25,6.2443e-25,6.2239e-25,6.2034e-25,6.1829e-25,6.1622e-25,6.1414e-25,6.1205e-25,6.0996e-25,6.0785e-25,6.0574e-25,6.0361e-25,6.0148e-25,5.9934e-25,5.9719e-25,5.9503e-25,5.9287e-25,5.9069e-25,5.8847e-25},
      std::array<float, 1001>{ 1.4466e-21,1.4445e-21,1.4424e-21,1.4403e-21,1.4382e-21,1.4361e-21,1.434e-21,1.4319e-21,1.4298e-21,1.4277e-21,1.4256e-21,1.4235e-21,1.4214e-21,1.4194e-21,1.4173e-21,1.4152e-21,1.4131e-21,1.4111e-21,1.409e-21,1.4069e-21,1.4048e-21,1.4027e-21,1.4007e-21,1.3986e-21,1.3965e-21,1.3944e-21,1.3923e-21,1.3902e-21,1.3881e-21,1.386e-21,1.3839e-21,1.3818e-21,1.3797e-21,1.3776e-21,1.3755e-21,1.3733e-21,1.3712e-21,1.3691e-21,1.3669e-21,1.3648e-21,1.3626e-21,1.3605e-21,1.3583e-21,1.3561e-21,1.3539e-21,1.3517e-21,1.3495e-21,1.3473e-21,1.3451e-21,1.3429e-21,1.3407e-21,1.3384e-21,1.3362e-21,1.334e-21,1.3317e-21,1.3294e-21,1.3272e-21,1.3249e-21,1.3226e-21,1.3203e-21,1.318e-21,1.3156e-21,1.3133e-21,1.311e-21,1.3086e-21,1.3063e-21,1.3039e-21,1.3015e-21,1.2991e-21,1.2968e-21,1.2944e-21,1.2919e-21,1.2895e-21,1.2871e-21,1.2847e-21,1.2822e-21,1.2798e-21,1.2773e-21,1.2748e-21,1.2723e-21,1.2699e-21,1.2674e-21,1.2649e-21,1.2623e-21,1.2598e-21,1.2573e-21,1.2548e-21,1.2522e-21,1.2497e-21,1.2471e-21,1.2445e-21,1.242e-21,1.2394e-21,1.2368e-21,1.2342e-21,1.2316e-21,1.229e-21,1.2264e-21,1.2237e-21,1.2211e-21,1.2185e-21,1.2158e-21,1.2132e-21,1.2105e-21,1.2079e-21,1.2052e-21,1.2025e-21,1.1999e-21,1.1972e-21,1.1945e-21,1.1918e-21,1.1891e-21,1.1864e-21,1.1837e-21,1.181e-21,1.1783e-21,1.1756e-21,1.1729e-21,1.1701e-21,1.1674e-21,1.1647e-21,1.1619e-21,1.1592e-21,1.1565e-21,1.1537e-21,1.151e-21,1.1482e-21,1.1455e-21,1.1427e-21,1.1399e-21,1.1372e-21,1.1344e-21,1.1316e-21,1.1289e-21,1.1261e-21,1.1233e-21,1.1206e-21,1.1178e-21,1.115e-21,1.1122e-21,1.1095e-21,1.1067e-21,1.1039e-21,1.1011e-21,1.0983e-21,1.0955e-21,1.0928e-21,1.09e-21,1.0872e-21,1.0844e-21,1.0816e-21,1.0788e-21,1.076e-21,1.0733e-21,1.0705e-21,1.0677e-21,1.0649e-21,1.0621e-21,1.0593e-21,1.0565e-21,1.0537e-21,1.051e-21,1.0482e-21,1.0454e-21,1.0426e-21,1.0398e-21,1.037e-21,1.0343e-21,1.0315e-21,1.0287e-21,1.0259e-21,1.0231e-21,1.0204e-21,1.0176e-21,1.0148e-21,1.012e-21,1.0093e-21,1.0065e-21,1.0038e-21,1.0011e-21,9.9839e-22,9.9567e-22,9.9296e-22,9.9024e-22,9.8753e-22,9.8482e-22,9.8212e-22,9.7941e-22,9.7671e-22,9.74e-22,9.713e-22,9.686e-22,9.6591e-22,9.6321e-22,9.6052e-22,9.5783e-22,9.5514e-22,9.5246e-22,9.4977e-22,9.4709e-22,9.4441e-22,9.4173e-22,9.3906e-22,9.3639e-22,9.3372e-22,9.3105e-22,9.2839e-22,9.2573e-22,9.2307e-22,9.2041e-22,9.1776e-22,9.1511e-22,9.1246e-22,9.0981e-22,9.0717e-22,9.0453e-22,9.0189e-22,8.9926e-22,8.9663e-22,8.94e-22,8.9138e-22,8.8875e-22,8.8614e-22,8.8352e-22,8.8091e-22,8.783e-22,8.7569e-22,8.7309e-22,8.7049e-22,8.679e-22,8.6531e-22,8.6272e-22,8.6013e-22,8.5755e-22,8.5498e-22,8.524e-22,8.4983e-22,8.4726e-22,8.447e-22,8.4214e-22,8.3959e-22,8.3703e-22,8.3449e-22,8.3194e-22,8.294e-22,8.2687e-22,8.2433e-22,8.2181e-22,8.1928e-22,8.1676e-22,8.1425e-22,8.1173e-22,8.0923e-22,8.0672e-22,8.0422e-22,8.0173e-22,7.9924e-22,7.9675e-22,7.9427e-22,7.9179e-22,7.8932e-22,7.8685e-22,7.8439e-22,7.8193e-22,7.7947e-22,7.7702e-22,7.7457e-22,7.7213e-22,7.6969e-22,7.6726e-22,7.6483e-22,7.6241e-22,7.5999e-22,7.5758e-22,7.5517e-22,7.5276e-22,7.5036e-22,7.4797e-22,7.4558e-22,7.432e-22,7.4082e-22,7.3844e-22,7.3607e-22,7.3371e-22,7.3135e-22,7.2899e-22,7.2664e-22,7.243e-22,7.2196e-22,7.1962e-22,7.1729e-22,7.1497e-22,7.1265e-22,7.1033e-22,7.0802e-22,7.0572e-22,7.0342e-22,7.0113e-22,6.9884e-22,6.9656e-22,6.9428e-22,6.9201e-22,6.8979e-22,6.8757e-22,6.8535e-22,6.8315e-22,6.8095e-22,6.7875e-22,6.7656e-22,6.7437e-22,6.7219e-22,6.7001e-22,6.6784e-22,6.6568e-22,6.6352e-22,6.6136e-22,6.5922e-22,6.5707e-22,6.5493e-22,6.528e-22,6.5067e-22,6.4855e-22,6.4643e-22,6.4432e-22,6.4222e-22,6.4012e-22,6.3802e-22,6.3593e-22,6.3385e-22,6.3177e-22,6.297e-22,6.2763e-22,6.2557e-22,6.2351e-22,6.2146e-22,6.1942e-22,6.1738e-22,6.1534e-22,6.1331e-22,6.1129e-22,6.0927e-22,6.0726e-22,6.0526e-22,6.0326e-22,6.0126e-22,5.9927e-22,5.9729e-22,5.9531e-22,5.9333e-22,5.9137e-22,5.8941e-22,5.8745e-22,5.855e-22,5.8355e-22,5.8161e-22,5.7968e-22,5.7775e-22,5.7583e-22,5.7391e-22,5.72e-22,5.701e-22,5.682e-22,5.663e-22,5.6441e-22,5.6253e-22,5.6065e-22,5.5878e-22,5.5691e-22,5.5505e-22,5.5319e-22,5.5134e-22,5.495e-22,5.4766e-22,5.4583e-22,5.44e-22,5.4218e-22,5.4036e-22,5.3855e-22,5.3674e-22,5.3494e-22,5.3315e-22,5.3136e-22,5.2958e-22,5.278e-22,5.2602e-22,5.2426e-22,5.225e-22,5.2074e-22,5.1899e-22,5.1724e-22,5.155e-22,5.1377e-22,5.1204e-22,5.1032e-22,5.086e-22,5.0689e-22,5.0518e-22,5.0348e-22,5.0178e-22,5.0013e-22,4.9849e-22,4.9685e-22,4.9521e-22,4.9358e-22,4.9196e-22,4.9034e-22,4.8872e-22,4.8712e-22,4.8551e-22,4.8391e-22,4.8232e-22,4.8073e-22,4.7915e-22,4.7757e-22,4.76e-22,4.7443e-22,4.7287e-22,4.7131e-22,4.6975e-22,4.6821e-22,4.6666e-22,4.6513e-22,4.6359e-22,4.6207e-22,4.6054e-22,4.5902e-22,4.5751e-22,4.56e-22,4.545e-22,4.53e-22,4.5151e-22,4.5002e-22,4.4854e-22,4.4706e-22,4.4559e-22,4.4412e-22,4.4265e-22,4.4119e-22,4.3974e-22,4.3829e-22,4.3684e-22,4.354e-22,4.3397e-22,4.3254e-22,4.3111e-22,4.2969e-22,4.2828e-22,4.2686e-22,4.2546e-22,4.2406e-22,4.2266e-22,4.2126e-22,4.1988e-22,4.1849e-22,4.1711e-22,4.1574e-22,4.1437e-22,4.13e-22,4.1164e-22,4.1029e-22,4.0894e-22,4.0759e-22,4.0625e-22,4.0491e-22,4.0357e-22,4.0225e-22,4.0092e-22,3.996e-22,3.9828e-22,3.9697e-22,3.9567e-22,3.9436e-22,3.9306e-22,3.9177e-22,3.9048e-22,3.892e-22,3.8791e-22,3.8664e-22,3.854e-22,3.8417e-22,3.8294e-22,3.8172e-22,3.8051e-22,3.7929e-22,3.7808e-22,3.7688e-22,3.7568e-22,3.7448e-22,3.7329e-22,3.721e-22,3.7091e-22,3.6973e-22,3.6856e-22,3.6738e-22,3.6621e-22,3.6505e-22,3.6389e-22,3.6273e-22,3.6158e-22,3.6043e-22,3.5928e-22,3.5814e-22,3.57e-22,3.5587e-22,3.5474e-22,3.5361e-22,3.5249e-22,3.5137e-22,3.5025e-22,3.4914e-22,3.4803e-22,3.4693e-22,3.4583e-22,3.4473e-22,3.4364e-22,3.4255e-22,3.4147e-22,3.4038e-22,3.393e-22,3.3823e-22,3.3716e-22,3.3609e-22,3.3503e-22,3.3397e-22,3.3291e-22,3.3186e-22,3.3081e-22,3.2976e-22,3.2872e-22,3.2768e-22,3.2664e-22,3.2561e-22,3.2458e-22,3.2356e-22,3.2253e-22,3.2151e-22,3.205e-22,3.1949e-22,3.1848e-22,3.1747e-22,3.1647e-22,3.1547e-22,3.1448e-22,3.1349e-22,3.125e-22,3.1155e-22,3.106e-22,3.0966e-22,3.0872e-22,3.0779e-22,3.0686e-22,3.0593e-22,3.05e-22,3.0408e-22,3.0316e-22,3.0224e-22,3.0133e-22,3.0042e-22,2.9951e-22,2.986e-22,2.977e-22,2.968e-22,2.9591e-22,2.9501e-22,2.9412e-22,2.9324e-22,2.9235e-22,2.9147e-22,2.9059e-22,2.8972e-22,2.8885e-22,2.8798e-22,2.8711e-22,2.8625e-22,2.8539e-22,2.8453e-22,2.8367e-22,2.8282e-22,2.8197e-22,2.8112e-22,2.8028e-22,2.7944e-22,2.786e-22,2.7776e-22,2.7693e-22,2.761e-22,2.7527e-22,2.7445e-22,2.7362e-22,2.728e-22,2.7199e-22,2.7117e-22,2.7036e-22,2.6955e-22,2.6875e-22,2.6794e-22,2.6714e-22,2.6634e-22,2.6555e-22,2.6475e-22,2.6396e-22,2.6317e-22,2.6239e-22,2.6164e-22,2.609e-22,2.6015e-22,2.5941e-22,2.5868e-22,2.5794e-22,2.5721e-22,2.5648e-22,2.5575e-22,2.5503e-22,2.5431e-22,2.5359e-22,2.5287e-22,2.5215e-22,2.5144e-22,2.5073e-22,2.5002e-22,2.4931e-22,2.4861e-22,2.479e-22,2.472e-22,2.4651e-22,2.4581e-22,2.4512e-22,2.4443e-22,2.4374e-22,2.4305e-22,2.4237e-22,2.4169e-22,2.4101e-22,2.4033e-22,2.3965e-22,2.3898e-22,2.3831e-22,2.3764e-22,2.3697e-22,2.3631e-22,2.3564e-22,2.3498e-22,2.3432e-22,2.3367e-22,2.3301e-22,2.3236e-22,2.3171e-22,2.3106e-22,2.3041e-22,2.2977e-22,2.2913e-22,2.2849e-22,2.2785e-22,2.2721e-22,2.2668e-22,2.2617e-22,2.2567e-22,2.2516e-22,2.2466e-22,2.2416e-22,2.2366e-22,2.2317e-22,2.2267e-22,2.2218e-22,2.2169e-22,2.2119e-22,2.2071e-22,2.2022e-22,2.1973e-22,2.1925e-22,2.1876e-22,2.1828e-22,2.178e-22,2.1733e-22,2.1685e-22,2.1637e-22,2.159e-22,2.1543e-22,2.1496e-22,2.1449e-22,2.1402e-22,2.1355e-22,2.1309e-22,2.1262e-22,2.1216e-22,2.117e-22,2.1124e-22,2.1078e-22,2.1033e-22,2.0987e-22,2.0942e-22,2.0897e-22,2.0852e-22,2.0807e-22,2.0762e-22,2.0717e-22,2.0673e-22,2.0629e-22,2.0584e-22,2.0539e-22,2.0475e-22,2.041e-22,2.0346e-22,2.0282e-22,2.0218e-22,2.0154e-22,2.009e-22,2.0027e-22,1.9963e-22,1.99e-22,1.9838e-22,1.9775e-22,1.9712e-22,1.965e-22,1.9588e-22,1.9526e-22,1.9464e-22,1.9403e-22,1.9341e-22,1.928e-22,1.9219e-22,1.9158e-22,1.9098e-22,1.9037e-22,1.8977e-22,1.8917e-22,1.8857e-22,1.8797e-22,1.8737e-22,1.8678e-22,1.8619e-22,1.856e-22,1.8501e-22,1.8442e-22,1.8384e-22,1.8325e-22,1.8267e-22,1.8209e-22,1.8151e-22,1.8093e-22,1.8036e-22,1.799e-22,1.7949e-22,1.7909e-22,1.787e-22,1.783e-22,1.779e-22,1.7751e-22,1.7711e-22,1.7672e-22,1.7633e-22,1.7594e-22,1.7555e-22,1.7516e-22,1.7477e-22,1.7439e-22,1.74e-22,1.7362e-22,1.7323e-22,1.7285e-22,1.7247e-22,1.7209e-22,1.7172e-22,1.7134e-22,1.7096e-22,1.7059e-22,1.7021e-22,1.6984e-22,1.6947e-22,1.691e-22,1.6873e-22,1.6836e-22,1.68e-22,1.6763e-22,1.6726e-22,1.669e-22,1.6654e-22,1.6618e-22,1.6582e-22,1.6548e-22,1.6515e-22,1.6482e-22,1.645e-22,1.6417e-22,1.6385e-22,1.6352e-22,1.632e-22,1.6288e-22,1.6256e-22,1.6224e-22,1.6192e-22,1.616e-22,1.6128e-22,1.6096e-22,1.6065e-22,1.6033e-22,1.6002e-22,1.5971e-22,1.5939e-22,1.5908e-22,1.5877e-22,1.5846e-22,1.5816e-22,1.5785e-22,1.5754e-22,1.5724e-22,1.5693e-22,1.5663e-22,1.5632e-22,1.5602e-22,1.5572e-22,1.5542e-22,1.5512e-22,1.5481e-22,1.5445e-22,1.5409e-22,1.5373e-22,1.5337e-22,1.5301e-22,1.5265e-22,1.523e-22,1.5194e-22,1.5159e-22,1.5124e-22,1.5088e-22,1.5053e-22,1.5018e-22,1.4983e-22,1.4948e-22,1.4914e-22,1.4879e-22,1.4844e-22,1.481e-22,1.4776e-22,1.4741e-22,1.4707e-22,1.4673e-22,1.4639e-22,1.4605e-22,1.4572e-22,1.4538e-22,1.4504e-22,1.4471e-22,1.4437e-22,1.4404e-22,1.4371e-22,1.4346e-22,1.4323e-22,1.43e-22,1.4277e-22,1.4253e-22,1.423e-22,1.4207e-22,1.4184e-22,1.4161e-22,1.4139e-22,1.4116e-22,1.4093e-22,1.407e-22,1.4048e-22,1.4025e-22,1.4003e-22,1.398e-22,1.3958e-22,1.3936e-22,1.3914e-22,1.3892e-22,1.3869e-22,1.3847e-22,1.3825e-22,1.3804e-22,1.3782e-22,1.376e-22,1.3738e-22,1.3717e-22,1.3695e-22,1.3722e-22,1.3751e-22,1.378e-22,1.381e-22,1.3839e-22,1.3868e-22,1.3897e-22,1.3926e-22,1.3955e-22,1.3984e-22,1.4013e-22,1.4042e-22,1.4071e-22,1.41e-22,1.4129e-22,1.4158e-22,1.4186e-22,1.4215e-22,1.4244e-22,1.4272e-22,1.4301e-22,1.433e-22,1.4358e-22,1.4387e-22,1.4415e-22,1.4444e-22,1.4472e-22,1.4501e-22,1.4469e-22,1.4431e-22,1.4394e-22,1.4358e-22,1.4321e-22,1.4284e-22,1.4248e-22,1.4211e-22,1.4175e-22,1.4138e-22,1.4102e-22,1.4066e-22,1.403e-22,1.3994e-22,1.3958e-22,1.3923e-22,1.3887e-22,1.3851e-22,1.3816e-22,1.3781e-22,1.3745e-22,1.371e-22,1.3675e-22,1.364e-22,1.3605e-22,1.357e-22,1.3536e-22,1.3501e-22,1.3467e-22,1.3432e-22,1.3398e-22,1.3363e-22,1.3329e-22,1.3295e-22,1.3261e-22,1.3227e-22,1.3193e-22,1.316e-22,1.3126e-22,1.3092e-22,1.3059e-22,1.3026e-22,1.2992e-22,1.2959e-22,1.2926e-22,1.2893e-22,1.286e-22,1.2827e-22,1.2794e-22,1.2761e-22,1.2729e-22,1.2674e-22,1.2613e-22,1.2552e-22,1.2492e-22,1.2432e-22,1.2373e-22,1.2314e-22,1.2256e-22,1.2198e-22,1.214e-22,1.2083e-22,1.2027e-22,1.1971e-22,1.1915e-22,1.186e-22,1.1806e-22,1.1752e-22,1.1698e-22,1.1645e-22,1.1592e-22,1.154e-22,1.1488e-22,1.1437e-22,1.1386e-22,1.1335e-22,1.1286e-22,1.1236e-22,1.1187e-22,1.1139e-22,1.1091e-22,1.1043e-22,1.0996e-22,1.0949e-22,1.0903e-22,1.0857e-22,1.0812e-22,1.0767e-22,1.0723e-22,1.0679e-22,1.0635e-22,1.0592e-22,1.055e-22,1.0508e-22,1.0466e-22,1.0425e-22,1.0385e-22} };

    const static float base_density = 1.3681e+20f; //km^-3

    const static float Hp = 730e-3f; // Continental value

    const static float pb = 2e6f; // background value in km^-3

    template <typename Float, typename UInt32, typename Mask, typename Spectrum, typename Wavelength>
    class UrbanAerosol : public GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength> {
        const float pb_bd;

    public:
        UrbanAerosol() : GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength>(tabulatedValues),
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
#endif //_AEROSOL_URBAN