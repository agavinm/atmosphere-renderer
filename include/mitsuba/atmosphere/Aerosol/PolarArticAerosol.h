#ifndef _AEROSOL_POLAR_ARTIC
#define _AEROSOL_POLAR_ARTIC

#include <array>
#include <mitsuba/atmosphere/GlobalAtmosphericAerosol.h>

namespace polarArticAerosol {

	// row[0] = wavelength
	// row[1] = cross section absorption coefficient
	// row[2] = cross section scattering coefficient
    const static std::array<std::array<float, 1001>, 3> tabulatedValues =
	{ std::array<float, 1001>{ 100,100.23,100.46,100.69,100.93,101.16,101.39,101.62,101.86,102.09,102.33,102.57,102.8,103.04,103.28,103.51,103.75,103.99,104.23,104.47,104.71,104.95,105.2,105.44,105.68,105.93,106.17,106.41,106.66,106.91,107.15,107.4,107.65,107.89,108.14,108.39,108.64,108.89,109.14,109.4,109.65,109.9,110.15,110.41,110.66,110.92,111.17,111.43,111.69,111.94,112.2,112.46,112.72,112.98,113.24,113.5,113.76,114.03,114.29,114.55,114.82,115.08,115.35,115.61,115.88,116.14,116.41,116.68,116.95,117.22,117.49,117.76,118.03,118.3,118.58,118.85,119.12,119.4,119.67,119.95,120.23,120.5,120.78,121.06,121.34,121.62,121.9,122.18,122.46,122.74,123.03,123.31,123.59,123.88,124.17,124.45,124.74,125.03,125.31,125.6,125.89,126.18,126.47,126.77,127.06,127.35,127.64,127.94,128.23,128.53,128.82,129.12,129.42,129.72,130.02,130.32,130.62,130.92,131.22,131.52,131.83,132.13,132.43,132.74,133.05,133.35,133.66,133.97,134.28,134.59,134.9,135.21,135.52,135.83,136.14,136.46,136.77,137.09,137.4,137.72,138.04,138.36,138.68,139,139.32,139.64,139.96,140.28,140.6,140.93,141.25,141.58,141.91,142.23,142.56,142.89,143.22,143.55,143.88,144.21,144.54,144.88,145.21,145.55,145.88,146.22,146.55,146.89,147.23,147.57,147.91,148.25,148.59,148.94,149.28,149.62,149.97,150.31,150.66,151.01,151.36,151.71,152.05,152.41,152.76,153.11,153.46,153.82,154.17,154.53,154.88,155.24,155.6,155.96,156.31,156.68,157.04,157.4,157.76,158.12,158.49,158.85,159.22,159.59,159.96,160.32,160.69,161.06,161.44,161.81,162.18,162.55,162.93,163.31,163.68,164.06,164.44,164.82,165.2,165.58,165.96,166.34,166.72,167.11,167.49,167.88,168.27,168.66,169.04,169.43,169.82,170.22,170.61,171,171.4,171.79,172.19,172.58,172.98,173.38,173.78,174.18,174.58,174.98,175.39,175.79,176.2,176.6,177.01,177.42,177.83,178.24,178.65,179.06,179.47,179.89,180.3,180.72,181.13,181.55,181.97,182.39,182.81,183.23,183.65,184.08,184.5,184.93,185.35,185.78,186.21,186.64,187.07,187.5,187.93,188.36,188.8,189.23,189.67,190.11,190.55,190.99,191.43,191.87,192.31,192.75,193.2,193.64,194.09,194.54,194.98,195.43,195.88,196.34,196.79,197.24,197.7,198.15,198.61,199.07,199.53,199.99,200.45,200.91,201.37,201.84,202.3,202.77,203.24,203.7,204.17,204.64,205.12,205.59,206.06,206.54,207.01,207.49,207.97,208.45,208.93,209.41,209.89,210.38,210.86,211.35,211.84,212.32,212.81,213.3,213.8,214.29,214.78,215.28,215.77,216.27,216.77,217.27,217.77,218.27,218.78,219.28,219.79,220.29,220.8,221.31,221.82,222.33,222.84,223.36,223.87,224.39,224.91,225.42,225.94,226.46,226.99,227.51,228.03,228.56,229.09,229.61,230.14,230.67,231.21,231.74,232.27,232.81,233.35,233.88,234.42,234.96,235.5,236.05,236.59,237.14,237.68,238.23,238.78,239.33,239.88,240.44,240.99,241.55,242.1,242.66,243.22,243.78,244.34,244.91,245.47,246.04,246.6,247.17,247.74,248.31,248.89,249.46,250.03,250.61,251.19,251.77,252.35,252.93,253.51,254.1,254.68,255.27,255.86,256.45,257.04,257.63,258.23,258.82,259.42,260.02,260.62,261.22,261.82,262.42,263.03,263.63,264.24,264.85,265.46,266.07,266.69,267.3,267.92,268.53,269.15,269.77,270.4,271.02,271.64,272.27,272.9,273.53,274.16,274.79,275.42,276.06,276.69,277.33,277.97,278.61,279.25,279.9,280.54,281.19,281.84,282.49,283.14,283.79,284.45,285.1,285.76,286.42,287.08,287.74,288.4,289.07,289.73,290.4,291.07,291.74,292.42,293.09,293.76,294.44,295.12,295.8,296.48,297.17,297.85,298.54,299.23,299.92,300.61,301.3,302,302.69,303.39,304.09,304.79,305.49,306.2,306.9,307.61,308.32,309.03,309.74,310.46,311.17,311.89,312.61,313.33,314.05,314.77,315.5,316.23,316.96,317.69,318.42,319.15,319.89,320.63,321.37,322.11,322.85,323.59,324.34,325.09,325.84,326.59,327.34,328.1,328.85,329.61,330.37,331.13,331.89,332.66,333.43,334.19,334.97,335.74,336.51,337.29,338.06,338.84,339.63,340.41,341.19,341.98,342.77,343.56,344.35,345.14,345.94,346.74,347.54,348.34,349.14,349.95,350.75,351.56,352.37,353.18,354,354.81,355.63,356.45,357.27,358.1,358.92,359.75,360.58,361.41,362.24,363.08,363.92,364.75,365.59,366.44,367.28,368.13,368.98,369.83,370.68,371.54,372.39,373.25,374.11,374.97,375.84,376.7,377.57,378.44,379.32,380.19,381.07,381.94,382.82,383.71,384.59,385.48,386.37,387.26,388.15,389.05,389.94,390.84,391.74,392.64,393.55,394.46,395.37,396.28,397.19,398.11,399.02,399.94,400.87,401.79,402.72,403.65,404.58,405.51,406.44,407.38,408.32,409.26,410.2,411.15,412.1,413.05,414,414.95,415.91,416.87,417.83,418.79,419.76,420.73,421.7,422.67,423.64,424.62,425.6,426.58,427.56,428.55,429.54,430.53,431.52,432.51,433.51,434.51,435.51,436.52,437.52,438.53,439.54,440.55,441.57,442.59,443.61,444.63,445.66,446.68,447.71,448.75,449.78,450.82,451.86,452.9,453.94,454.99,456.04,457.09,458.14,459.2,460.26,461.32,462.38,463.45,464.52,465.59,466.66,467.74,468.81,469.89,470.98,472.06,473.15,474.24,475.34,476.43,477.53,478.63,479.73,480.84,481.95,483.06,484.17,485.29,486.41,487.53,488.65,489.78,490.91,492.04,493.17,494.31,495.45,496.59,497.74,498.88,500.03,501.19,502.34,503.5,504.66,505.82,506.99,508.16,509.33,510.51,511.68,512.86,514.04,515.23,516.42,517.61,518.8,520,521.19,522.4,523.6,524.81,526.02,527.23,528.45,529.66,530.88,532.11,533.33,534.56,535.8,537.03,538.27,539.51,540.75,542,543.25,544.5,545.76,547.02,548.28,549.54,550.81,552.08,553.35,554.63,555.9,557.19,558.47,559.76,561.05,562.34,563.64,564.94,566.24,567.54,568.85,570.16,571.48,572.8,574.12,575.44,576.77,578.1,579.43,580.76,582.1,583.45,584.79,586.14,587.49,588.84,590.2,591.56,592.93,594.29,595.66,597.04,598.41,599.79,601.17,602.56,603.95,605.34,606.74,608.13,609.54,610.94,612.35,613.76,615.18,616.6,618.02,619.44,620.87,622.3,623.73,625.17,626.61,628.06,629.51,630.96,632.41,633.87,635.33,636.8,638.26,639.73,641.21,642.69,644.17,645.65,647.14,648.63,650.13,651.63,653.13,654.64,656.15,657.66,659.17,660.69,662.22,663.74,665.27,666.81,668.34,669.88,671.43,672.98,674.53,676.08,677.64,679.2,680.77,682.34,683.91,685.49,687.07,688.65,690.24,691.83,693.43,695.02,696.63,698.23,699.84,701.46,703.07,704.69,706.32,707.95,709.58,711.21,712.85,714.5,716.14,717.79,719.45,721.11,722.77,724.44,726.11,727.78,729.46,731.14,732.82,734.51,736.21,737.9,739.61,741.31,743.02,744.73,746.45,748.17,749.89,751.62,753.36,755.09,756.83,758.58,760.33,762.08,763.84,765.6,767.36,769.13,770.9,772.68,774.46,776.25,778.04,779.83,781.63,783.43,785.24,787.05,788.86,790.68,792.5,794.33,796.16,797.99,799.83,801.68,803.53,805.38,807.24,809.1,810.96,812.83,814.7,816.58,818.46,820.35,822.24,824.14,826.04,827.94,829.85,831.76,833.68,835.6,837.53,839.46,841.4,843.33,845.28,847.23,849.18,851.14,853.1,855.07,857.04,859.01,860.99,862.98,864.97,866.96,868.96,870.96,872.97,874.98,877,879.02,881.05,883.08,885.12,887.16,889.2,891.25,893.31,895.36,897.43,899.5,901.57,903.65,905.73,907.82,909.91,912.01,914.11,916.22,918.33,920.45,922.57,924.7,926.83,928.97,931.11,933.25,935.41,937.56,939.72,941.89,944.06,946.24,948.42,950.6,952.8,954.99,957.19,959.4,961.61,963.83,966.05,968.28,970.51,972.75,974.99,977.24,979.49,981.75,984.01,986.28,988.55,990.83,993.12,995.41,997.7,1000 },
      std::array<float, 1001>{ 9.3252e-17,9.3439e-17,9.3594e-17,9.3661e-17,9.3653e-17,9.3615e-17,9.3565e-17,9.3542e-17,9.3569e-17,9.3611e-17,9.3634e-17,9.359e-17,9.3495e-17,9.3326e-17,9.3162e-17,9.3003e-17,9.2946e-17,9.3025e-17,9.3186e-17,9.3453e-17,9.3712e-17,9.3921e-17,9.4026e-17,9.4003e-17,9.3879e-17,9.3664e-17,9.3454e-17,9.3296e-17,9.3209e-17,9.3217e-17,9.3298e-17,9.3333e-17,9.3385e-17,9.3372e-17,9.3346e-17,9.3287e-17,9.3292e-17,9.3347e-17,9.3483e-17,9.3664e-17,9.3846e-17,9.398e-17,9.4052e-17,9.3977e-17,9.3803e-17,9.3569e-17,9.3286e-17,9.3045e-17,9.2916e-17,9.2896e-17,9.298e-17,9.3186e-17,9.3408e-17,9.3596e-17,9.3757e-17,9.3807e-17,9.3785e-17,9.373e-17,9.3655e-17,9.3609e-17,9.3604e-17,9.3645e-17,9.3675e-17,9.3701e-17,9.3674e-17,9.3573e-17,9.342e-17,9.3216e-17,9.3042e-17,9.2929e-17,9.2922e-17,9.3026e-17,9.3212e-17,9.3519e-17,9.3773e-17,9.401e-17,9.4173e-17,9.4179e-17,9.4095e-17,9.3915e-17,9.3687e-17,9.3474e-17,9.3318e-17,9.3231e-17,9.3217e-17,9.3256e-17,9.3329e-17,9.34e-17,9.346e-17,9.3467e-17,9.3446e-17,9.3397e-17,9.3374e-17,9.3395e-17,9.349e-17,9.3625e-17,9.3817e-17,9.4022e-17,9.4171e-17,9.4249e-17,9.4242e-17,9.4097e-17,9.3854e-17,9.3594e-17,9.3311e-17,9.3075e-17,9.2944e-17,9.2918e-17,9.3042e-17,9.3205e-17,9.3452e-17,9.3692e-17,9.3892e-17,9.4017e-17,9.407e-17,9.4053e-17,9.3977e-17,9.392e-17,9.3865e-17,9.3826e-17,9.3843e-17,9.3882e-17,9.3928e-17,9.3962e-17,9.3955e-17,9.3885e-17,9.3782e-17,9.362e-17,9.3432e-17,9.3271e-17,9.3142e-17,9.3095e-17,9.315e-17,9.3315e-17,9.3577e-17,9.3876e-17,9.4156e-17,9.44e-17,9.4567e-17,9.4629e-17,9.4587e-17,9.4441e-17,9.4237e-17,9.3993e-17,9.3754e-17,9.3605e-17,9.348e-17,9.344e-17,9.3483e-17,9.3574e-17,9.367e-17,9.3751e-17,9.3792e-17,9.3821e-17,9.3775e-17,9.3732e-17,9.3705e-17,9.3676e-17,9.3712e-17,9.3797e-17,9.3948e-17,9.4142e-17,9.4348e-17,9.4546e-17,9.4675e-17,9.4726e-17,9.4681e-17,9.4551e-17,9.4325e-17,9.4058e-17,9.3753e-17,9.3503e-17,9.3281e-17,9.3143e-17,9.3098e-17,9.3174e-17,9.3333e-17,9.3562e-17,9.3842e-17,9.4122e-17,9.4345e-17,9.4504e-17,9.4576e-17,9.4589e-17,9.4513e-17,9.4383e-17,9.4248e-17,9.4137e-17,9.4049e-17,9.402e-17,9.4001e-17,9.4035e-17,9.4059e-17,9.409e-17,9.4055e-17,9.4009e-17,9.3905e-17,9.3791e-17,9.3647e-17,9.3548e-17,9.3463e-17,9.3492e-17,9.3562e-17,9.3728e-17,9.3968e-17,9.4241e-17,9.4519e-17,9.4778e-17,9.4958e-17,9.5039e-17,9.5036e-17,9.4893e-17,9.4682e-17,9.4403e-17,9.407e-17,9.3736e-17,9.3469e-17,9.3259e-17,9.3174e-17,9.3159e-17,9.3275e-17,9.348e-17,9.3717e-17,9.4e-17,9.4236e-17,9.4467e-17,9.4632e-17,9.4713e-17,9.473e-17,9.468e-17,9.4577e-17,9.4453e-17,9.4376e-17,9.427e-17,9.4234e-17,9.4212e-17,9.4227e-17,9.4247e-17,9.427e-17,9.4289e-17,9.428e-17,9.4233e-17,9.415e-17,9.4035e-17,9.3911e-17,9.3755e-17,9.3653e-17,9.3532e-17,9.3515e-17,9.3547e-17,9.3683e-17,9.3877e-17,9.4119e-17,9.4415e-17,9.4715e-17,9.4975e-17,9.5215e-17,9.5371e-17,9.5406e-17,9.5348e-17,9.518e-17,9.4934e-17,9.462e-17,9.4289e-17,9.3952e-17,9.3664e-17,9.3451e-17,9.3302e-17,9.3259e-17,9.333e-17,9.3455e-17,9.3682e-17,9.3935e-17,9.4199e-17,9.4444e-17,9.4657e-17,9.4808e-17,9.4907e-17,9.4967e-17,9.4948e-17,9.4891e-17,9.4808e-17,9.4711e-17,9.4638e-17,9.4541e-17,9.4531e-17,9.4505e-17,9.4518e-17,9.4548e-17,9.4598e-17,9.465e-17,9.4664e-17,9.4649e-17,9.4604e-17,9.4522e-17,9.4408e-17,9.4257e-17,9.4098e-17,9.3928e-17,9.3786e-17,9.3669e-17,9.364e-17,9.3661e-17,9.3755e-17,9.3951e-17,9.4211e-17,9.4528e-17,9.484e-17,9.5172e-17,9.5467e-17,9.5659e-17,9.5821e-17,9.5832e-17,9.5772e-17,9.5597e-17,9.5314e-17,9.5004e-17,9.465e-17,9.4291e-17,9.3953e-17,9.3663e-17,9.3467e-17,9.3385e-17,9.3348e-17,9.3438e-17,9.3646e-17,9.387e-17,9.4171e-17,9.4491e-17,9.4773e-17,9.5064e-17,9.5279e-17,9.5445e-17,9.554e-17,9.5551e-17,9.5522e-17,9.5415e-17,9.5251e-17,9.5123e-17,9.4941e-17,9.4791e-17,9.4688e-17,9.4545e-17,9.4506e-17,9.4472e-17,9.4476e-17,9.4509e-17,9.4577e-17,9.4583e-17,9.4627e-17,9.4659e-17,9.4625e-17,9.4611e-17,9.457e-17,9.4458e-17,9.4379e-17,9.4298e-17,9.4195e-17,9.4191e-17,9.416e-17,9.4221e-17,9.4325e-17,9.4488e-17,9.4674e-17,9.4913e-17,9.5194e-17,9.5408e-17,9.5635e-17,9.5855e-17,9.5953e-17,9.6051e-17,9.6072e-17,9.5929e-17,9.5785e-17,9.5538e-17,9.521e-17,9.4913e-17,9.4545e-17,9.416e-17,9.3886e-17,9.3599e-17,9.3362e-17,9.3305e-17,9.3261e-17,9.3325e-17,9.3523e-17,9.3727e-17,9.4021e-17,9.4382e-17,9.4751e-17,9.5089e-17,9.5448e-17,9.577e-17,9.6008e-17,9.6154e-17,9.6294e-17,9.6305e-17,9.6219e-17,9.6124e-17,9.591e-17,9.5613e-17,9.5391e-17,9.5086e-17,9.4759e-17,9.4545e-17,9.4328e-17,9.4074e-17,9.4009e-17,9.3955e-17,9.3895e-17,9.3987e-17,9.4113e-17,9.4205e-17,9.4396e-17,9.4587e-17,9.4729e-17,9.4878e-17,9.5053e-17,9.517e-17,9.5222e-17,9.5279e-17,9.5308e-17,9.5285e-17,9.5249e-17,9.5221e-17,9.5195e-17,9.5095e-17,9.5105e-17,9.5157e-17,9.5089e-17,9.5116e-17,9.5241e-17,9.5238e-17,9.5255e-17,9.5433e-17,9.5462e-17,9.5411e-17,9.5522e-17,9.5553e-17,9.5376e-17,9.5336e-17,9.5296e-17,9.5075e-17,9.486e-17,9.4778e-17,9.4548e-17,9.4311e-17,9.4179e-17,9.4059e-17,9.3936e-17,9.387e-17,9.3926e-17,9.3963e-17,9.4085e-17,9.4267e-17,9.4521e-17,9.4828e-17,9.5096e-17,9.5429e-17,9.5794e-17,9.61e-17,9.6315e-17,9.6608e-17,9.6866e-17,9.6862e-17,9.6876e-17,9.6973e-17,9.6851e-17,9.652e-17,9.6335e-17,9.6126e-17,9.5643e-17,9.522e-17,9.4997e-17,9.4596e-17,9.4063e-17,9.382e-17,9.3679e-17,9.3402e-17,9.3177e-17,9.3283e-17,9.337e-17,9.3391e-17,9.359e-17,9.3971e-17,9.4286e-17,9.4585e-17,9.5049e-17,9.5503e-17,9.5875e-17,9.6243e-17,9.6602e-17,9.6908e-17,9.7153e-17,9.7322e-17,9.7378e-17,9.739e-17,9.7405e-17,9.7233e-17,9.6953e-17,9.6757e-17,9.6571e-17,9.6128e-17,9.5687e-17,9.5482e-17,9.5252e-17,9.4788e-17,9.4432e-17,9.44e-17,9.4297e-17,9.3972e-17,9.3848e-17,9.4121e-17,9.4179e-17,9.4031e-17,9.422e-17,9.4669e-17,9.4838e-17,9.485e-17,9.5145e-17,9.5648e-17,9.5837e-17,9.5849e-17,9.6108e-17,9.6489e-17,9.6568e-17,9.6487e-17,9.659e-17,9.6797e-17,9.6737e-17,9.6568e-17,9.6546e-17,9.6549e-17,9.6437e-17,9.6285e-17,9.6215e-17,9.6121e-17,9.6013e-17,9.5977e-17,9.5957e-17,9.5854e-17,9.5784e-17,9.5877e-17,9.5989e-17,9.5891e-17,9.58e-17,9.5999e-17,9.6201e-17,9.6106e-17,9.5947e-17,9.6115e-17,9.638e-17,9.6296e-17,9.5995e-17,9.6008e-17,9.6283e-17,9.6288e-17,9.5854e-17,9.5654e-17,9.5885e-17,9.6006e-17,9.5622e-17,9.5229e-17,9.5395e-17,9.5663e-17,9.5509e-17,9.5106e-17,9.5074e-17,9.5481e-17,9.5673e-17,9.5433e-17,9.5328e-17,9.5724e-17,9.6218e-17,9.6305e-17,9.6199e-17,9.6481e-17,9.7015e-17,9.7364e-17,9.7391e-17,9.7503e-17,9.7893e-17,9.8298e-17,9.8362e-17,9.8347e-17,9.8483e-17,9.8671e-17,9.8738e-17,9.853e-17,9.8424e-17,9.8361e-17,9.8228e-17,9.7925e-17,9.7656e-17,9.746e-17,9.7146e-17,9.6685e-17,9.6279e-17,9.6091e-17,9.5907e-17,9.5508e-17,9.4971e-17,9.4744e-17,9.4823e-17,9.4826e-17,9.4479e-17,9.4133e-17,9.4224e-17,9.4661e-17,9.4895e-17,9.4755e-17,9.4588e-17,9.4969e-17,9.5765e-17,9.6303e-17,9.6292e-17,9.6269e-17,9.6836e-17,9.7802e-17,9.851e-17,9.8546e-17,9.8485e-17,9.8948e-17,9.9887e-17,1.0058e-16,1.0056e-16,1.0026e-16,1.0041e-16,1.011e-16,1.0163e-16,1.0147e-16,1.0089e-16,1.0067e-16,1.0092e-16,1.0122e-16,1.0095e-16,1.0016e-16,9.9553e-17,9.949e-17,9.9595e-17,9.9332e-17,9.8562e-17,9.7762e-17,9.7458e-17,9.7495e-17,9.7394e-17,9.686e-17,9.6199e-17,9.5858e-17,9.5925e-17,9.6038e-17,9.5943e-17,9.5599e-17,9.5449e-17,9.5591e-17,9.5882e-17,9.6114e-17,9.6149e-17,9.6174e-17,9.6415e-17,9.6807e-17,9.7259e-17,9.7624e-17,9.7836e-17,9.8108e-17,9.8527e-17,9.9048e-17,9.9492e-17,9.981e-17,1.0009e-16,1.0044e-16,1.0089e-16,1.0134e-16,1.0158e-16,1.0179e-16,1.0198e-16,1.0225e-16,1.0252e-16,1.0271e-16,1.0272e-16,1.0267e-16,1.0273e-16,1.0282e-16,1.0291e-16,1.0286e-16,1.0259e-16,1.0231e-16,1.0219e-16,1.0221e-16,1.0217e-16,1.0193e-16,1.0148e-16,1.0104e-16,1.0082e-16,1.0083e-16,1.0082e-16,1.0054e-16,1.0002e-16,9.9518e-17,9.9251e-17,9.9339e-17,9.948e-17,9.937e-17,9.896e-17,9.8455e-17,9.8181e-17,9.8361e-17,9.874e-17,9.8909e-17,9.8712e-17,9.8262e-17,9.7979e-17,9.8139e-17,9.8697e-17,9.9184e-17,9.9303e-17,9.9029e-17,9.8667e-17,9.873e-17,9.9262e-17,1.0004e-16,1.0049e-16,1.0048e-16,1.0016e-16,9.9906e-17,1.0025e-16,1.0104e-16,1.0183e-16,1.0221e-16,1.0205e-16,1.0163e-16,1.0154e-16,1.0206e-16,1.0297e-16,1.0369e-16,1.0392e-16,1.036e-16,1.0314e-16,1.0308e-16,1.0369e-16,1.0457e-16,1.0522e-16,1.053e-16,1.0484e-16,1.0423e-16,1.0411e-16,1.0468e-16,1.0551e-16,1.0609e-16,1.0607e-16,1.0549e-16,1.0474e-16,1.0449e-16,1.0493e-16,1.057e-16,1.0624e-16,1.0618e-16,1.0554e-16,1.0463e-16,1.0418e-16,1.0444e-16,1.0506e-16,1.0562e-16,1.0562e-16,1.05e-16,1.0404e-16,1.0336e-16,1.0333e-16,1.0385e-16,1.0441e-16,1.0463e-16,1.0417e-16,1.032e-16,1.0233e-16,1.0205e-16,1.024e-16,1.0306e-16,1.0352e-16,1.0342e-16,1.0265e-16,1.0165e-16,1.0105e-16,1.0114e-16,1.0178e-16,1.025e-16,1.0287e-16,1.0262e-16,1.0183e-16,1.0104e-16,1.0079e-16,1.0123e-16,1.0208e-16,1.0294e-16,1.0337e-16,1.0311e-16,1.0243e-16,1.0191e-16,1.0193e-16,1.026e-16,1.0364e-16,1.0461e-16,1.0509e-16,1.049e-16,1.0435e-16,1.0397e-16,1.0413e-16,1.0488e-16,1.0599e-16,1.0698e-16,1.0751e-16,1.073e-16,1.0676e-16,1.0636e-16,1.0641e-16,1.0707e-16,1.0809e-16,1.0905e-16,1.0955e-16,1.0936e-16,1.0869e-16,1.0812e-16,1.0791e-16,1.0825e-16,1.0898e-16,1.0981e-16,1.1031e-16,1.102e-16,1.0943e-16,1.086e-16,1.0805e-16,1.0793e-16,1.0828e-16,1.0886e-16,1.0939e-16,1.0942e-16,1.0881e-16,1.0784e-16,1.07e-16,1.0652e-16,1.0643e-16,1.0668e-16,1.0709e-16,1.0739e-16,1.0721e-16,1.0642e-16,1.0547e-16,1.0477e-16,1.044e-16,1.0431e-16,1.045e-16,1.0487e-16,1.0518e-16,1.0499e-16,1.0426e-16,1.0349e-16,1.0301e-16,1.0282e-16,1.0283e-16,1.031e-16,1.0354e-16,1.0405e-16,1.0411e-16,1.0363e-16,1.0309e-16,1.0285e-16,1.0284e-16,1.0295e-16,1.0324e-16,1.0381e-16,1.0454e-16,1.0498e-16,1.0483e-16,1.0447e-16,1.044e-16,1.045e-16,1.0467e-16,1.05e-16,1.0559e-16,1.065e-16,1.0732e-16,1.0769e-16,1.0764e-16,1.0763e-16,1.0793e-16,1.0828e-16,1.0848e-16,1.0879e-16,1.0947e-16,1.1038e-16,1.1105e-16,1.1126e-16,1.1125e-16,1.1139e-16,1.1176e-16,1.1198e-16,1.1197e-16,1.1211e-16,1.127e-16,1.1343e-16,1.1387e-16,1.1391e-16,1.1387e-16,1.1401e-16,1.1432e-16,1.142e-16,1.1375e-16,1.1345e-16,1.1364e-16,1.1399e-16,1.1399e-16,1.1371e-16,1.1345e-16,1.1333e-16,1.134e-16,1.1322e-16,1.1263e-16,1.1209e-16,1.1204e-16,1.1217e-16,1.1209e-16,1.1179e-16,1.1154e-16,1.1147e-16,1.1157e-16,1.116e-16,1.1119e-16,1.1054e-16,1.1022e-16,1.1024e-16,1.1016e-16,1.0995e-16,1.0978e-16,1.0982e-16,1.1005e-16,1.1041e-16,1.1051e-16,1.1015e-16,1.0964e-16,1.0942e-16,1.0939e-16,1.0931e-16,1.092e-16,1.0929e-16,1.0962e-16,1.1021e-16,1.1084e-16,1.1119e-16,1.1104e-16,1.1063e-16,1.1041e-16,1.1029e-16,1.1024e-16,1.1021e-16,1.1042e-16,1.1093e-16,1.1172e-16,1.1265e-16,1.1338e-16,1.1348e-16,1.1312e-16,1.1276e-16,1.1256e-16,1.1239e-16,1.1228e-16,1.1243e-16,1.1288e-16,1.1376e-16,1.1472e-16,1.1567e-16,1.1612e-16,1.1591e-16,1.154e-16,1.1492e-16,1.1453e-16,1.1414e-16,1.1395e-16,1.1407e-16,1.1457e-16,1.1547e-16,1.165e-16,1.1728e-16,1.1748e-16,1.1705e-16,1.1633e-16,1.1569e-16,1.1503e-16,1.1435e-16,1.1393e-16,1.1393e-16,1.1428e-16,1.1507e-16,1.1597e-16,1.167e-16,1.1683e-16,1.1636e-16,1.1561e-16,1.1479e-16,1.1389e-16,1.1298e-16,1.1232e-16,1.1204e-16},
      std::array<float, 1001>{ 2.6076e-17,2.6077e-17,2.6078e-17,2.6079e-17,2.6081e-17,2.6082e-17,2.6083e-17,2.6084e-17,2.6086e-17,2.6087e-17,2.6088e-17,2.609e-17,2.6091e-17,2.6093e-17,2.6094e-17,2.6095e-17,2.6097e-17,2.6098e-17,2.6099e-17,2.61e-17,2.6102e-17,2.6103e-17,2.6104e-17,2.6106e-17,2.6107e-17,2.6109e-17,2.611e-17,2.6112e-17,2.6113e-17,2.6114e-17,2.6116e-17,2.6117e-17,2.6118e-17,2.612e-17,2.6121e-17,2.6122e-17,2.6123e-17,2.6125e-17,2.6126e-17,2.6128e-17,2.6129e-17,2.6131e-17,2.6132e-17,2.6134e-17,2.6135e-17,2.6137e-17,2.6138e-17,2.6139e-17,2.6141e-17,2.6142e-17,2.6143e-17,2.6144e-17,2.6146e-17,2.6147e-17,2.6149e-17,2.615e-17,2.6152e-17,2.6153e-17,2.6155e-17,2.6156e-17,2.6158e-17,2.6159e-17,2.6161e-17,2.6162e-17,2.6164e-17,2.6165e-17,2.6166e-17,2.6168e-17,2.6169e-17,2.617e-17,2.6172e-17,2.6173e-17,2.6175e-17,2.6177e-17,2.6178e-17,2.618e-17,2.6182e-17,2.6183e-17,2.6185e-17,2.6186e-17,2.6188e-17,2.6189e-17,2.619e-17,2.6192e-17,2.6193e-17,2.6194e-17,2.6196e-17,2.6197e-17,2.6199e-17,2.6201e-17,2.6202e-17,2.6204e-17,2.6206e-17,2.6207e-17,2.6209e-17,2.6211e-17,2.6212e-17,2.6214e-17,2.6215e-17,2.6216e-17,2.6218e-17,2.6219e-17,2.6221e-17,2.6222e-17,2.6224e-17,2.6225e-17,2.6227e-17,2.6228e-17,2.623e-17,2.6232e-17,2.6234e-17,2.6235e-17,2.6237e-17,2.6239e-17,2.624e-17,2.6242e-17,2.6244e-17,2.6245e-17,2.6246e-17,2.6248e-17,2.6249e-17,2.6251e-17,2.6252e-17,2.6254e-17,2.6255e-17,2.6257e-17,2.6259e-17,2.626e-17,2.6262e-17,2.6264e-17,2.6266e-17,2.6268e-17,2.627e-17,2.6271e-17,2.6273e-17,2.6275e-17,2.6276e-17,2.6278e-17,2.6279e-17,2.628e-17,2.6282e-17,2.6283e-17,2.6285e-17,2.6288e-17,2.6288e-17,2.629e-17,2.6292e-17,2.6294e-17,2.6296e-17,2.6298e-17,2.63e-17,2.6301e-17,2.6303e-17,2.6305e-17,2.6307e-17,2.6308e-17,2.631e-17,2.6311e-17,2.6313e-17,2.6314e-17,2.6316e-17,2.6317e-17,2.6319e-17,2.6321e-17,2.6323e-17,2.6324e-17,2.6326e-17,2.6328e-17,2.633e-17,2.6332e-17,2.6334e-17,2.6336e-17,2.6338e-17,2.634e-17,2.6342e-17,2.6344e-17,2.6345e-17,2.6347e-17,2.6349e-17,2.635e-17,2.6351e-17,2.6353e-17,2.6355e-17,2.6356e-17,2.6358e-17,2.636e-17,2.6362e-17,2.6364e-17,2.6366e-17,2.6368e-17,2.637e-17,2.6373e-17,2.6375e-17,2.6377e-17,2.6379e-17,2.6381e-17,2.6382e-17,2.6384e-17,2.6386e-17,2.6387e-17,2.6389e-17,2.639e-17,2.6392e-17,2.6393e-17,2.6396e-17,2.6397e-17,2.6399e-17,2.6401e-17,2.6403e-17,2.6405e-17,2.6407e-17,2.6409e-17,2.6412e-17,2.6414e-17,2.6416e-17,2.6418e-17,2.642e-17,2.6423e-17,2.6425e-17,2.6427e-17,2.6428e-17,2.643e-17,2.6431e-17,2.6433e-17,2.6435e-17,2.6436e-17,2.6438e-17,2.644e-17,2.6442e-17,2.6444e-17,2.6446e-17,2.6448e-17,2.645e-17,2.6452e-17,2.6455e-17,2.6457e-17,2.6459e-17,2.6461e-17,2.6464e-17,2.6466e-17,2.6468e-17,2.6471e-17,2.6473e-17,2.6475e-17,2.6476e-17,2.6478e-17,2.648e-17,2.6481e-17,2.6483e-17,2.6485e-17,2.6486e-17,2.6488e-17,2.6491e-17,2.6493e-17,2.6495e-17,2.6497e-17,2.6499e-17,2.6502e-17,2.6504e-17,2.6506e-17,2.6509e-17,2.6511e-17,2.6514e-17,2.6516e-17,2.6519e-17,2.6521e-17,2.6523e-17,2.6525e-17,2.6527e-17,2.6529e-17,2.6531e-17,2.6532e-17,2.6534e-17,2.6536e-17,2.6538e-17,2.654e-17,2.6542e-17,2.6544e-17,2.6546e-17,2.6549e-17,2.6551e-17,2.6553e-17,2.6556e-17,2.6559e-17,2.656e-17,2.6563e-17,2.6565e-17,2.6568e-17,2.6571e-17,2.6573e-17,2.6576e-17,2.6579e-17,2.6581e-17,2.6583e-17,2.6585e-17,2.6587e-17,2.6589e-17,2.6591e-17,2.6593e-17,2.6594e-17,2.6596e-17,2.6598e-17,2.6601e-17,2.6603e-17,2.6605e-17,2.6608e-17,2.661e-17,2.6612e-17,2.6615e-17,2.6617e-17,2.6619e-17,2.6621e-17,2.6624e-17,2.6626e-17,2.6629e-17,2.6632e-17,2.6635e-17,2.6638e-17,2.6641e-17,2.6643e-17,2.6646e-17,2.6648e-17,2.6651e-17,2.6652e-17,2.6654e-17,2.6656e-17,2.6658e-17,2.666e-17,2.6663e-17,2.6665e-17,2.6667e-17,2.667e-17,2.6672e-17,2.6675e-17,2.6677e-17,2.668e-17,2.6682e-17,2.6684e-17,2.6686e-17,2.6688e-17,2.6691e-17,2.6693e-17,2.6695e-17,2.6698e-17,2.6701e-17,2.6704e-17,2.6707e-17,2.671e-17,2.6713e-17,2.6716e-17,2.6719e-17,2.6721e-17,2.6723e-17,2.6726e-17,2.6728e-17,2.673e-17,2.6732e-17,2.6734e-17,2.6737e-17,2.6739e-17,2.6742e-17,2.6745e-17,2.6748e-17,2.675e-17,2.6753e-17,2.6755e-17,2.6758e-17,2.676e-17,2.6762e-17,2.6764e-17,2.6766e-17,2.6768e-17,2.6771e-17,2.6773e-17,2.6776e-17,2.6779e-17,2.6782e-17,2.6785e-17,2.6788e-17,2.6791e-17,2.6794e-17,2.6797e-17,2.68e-17,2.6803e-17,2.6805e-17,2.6808e-17,2.681e-17,2.6813e-17,2.6815e-17,2.6818e-17,2.682e-17,2.6823e-17,2.6826e-17,2.6829e-17,2.6832e-17,2.6834e-17,2.6837e-17,2.684e-17,2.6842e-17,2.6845e-17,2.6847e-17,2.6849e-17,2.6851e-17,2.6853e-17,2.6855e-17,2.6857e-17,2.686e-17,2.6862e-17,2.6865e-17,2.6868e-17,2.6871e-17,2.6874e-17,2.6877e-17,2.688e-17,2.6884e-17,2.6887e-17,2.689e-17,2.6893e-17,2.6896e-17,2.6899e-17,2.6902e-17,2.6904e-17,2.6907e-17,2.691e-17,2.6913e-17,2.6915e-17,2.6918e-17,2.6921e-17,2.6924e-17,2.6927e-17,2.693e-17,2.6933e-17,2.6936e-17,2.6939e-17,2.6941e-17,2.6943e-17,2.6945e-17,2.6947e-17,2.695e-17,2.6951e-17,2.6953e-17,2.6955e-17,2.6957e-17,2.696e-17,2.6962e-17,2.6965e-17,2.6968e-17,2.6971e-17,2.6974e-17,2.6977e-17,2.698e-17,2.6983e-17,2.6986e-17,2.6989e-17,2.6992e-17,2.6995e-17,2.6998e-17,2.7001e-17,2.7005e-17,2.7007e-17,2.701e-17,2.7013e-17,2.7016e-17,2.7019e-17,2.7022e-17,2.7025e-17,2.7028e-17,2.7032e-17,2.7035e-17,2.7038e-17,2.704e-17,2.7043e-17,2.7046e-17,2.7049e-17,2.705e-17,2.7052e-17,2.7054e-17,2.7055e-17,2.7057e-17,2.7059e-17,2.706e-17,2.7062e-17,2.7064e-17,2.7066e-17,2.7068e-17,2.707e-17,2.7073e-17,2.7075e-17,2.7077e-17,2.708e-17,2.7083e-17,2.7085e-17,2.7087e-17,2.709e-17,2.7093e-17,2.7095e-17,2.7098e-17,2.71e-17,2.7103e-17,2.7105e-17,2.7108e-17,2.711e-17,2.7113e-17,2.7115e-17,2.7118e-17,2.712e-17,2.7122e-17,2.7125e-17,2.7127e-17,2.713e-17,2.7132e-17,2.7134e-17,2.7136e-17,2.7138e-17,2.714e-17,2.7142e-17,2.7143e-17,2.7144e-17,2.7145e-17,2.7146e-17,2.7147e-17,2.7147e-17,2.7148e-17,2.7149e-17,2.7149e-17,2.715e-17,2.7151e-17,2.7151e-17,2.7152e-17,2.7153e-17,2.7154e-17,2.7155e-17,2.7156e-17,2.7157e-17,2.7158e-17,2.7158e-17,2.7159e-17,2.716e-17,2.7161e-17,2.7162e-17,2.7163e-17,2.7164e-17,2.7164e-17,2.7165e-17,2.7166e-17,2.7166e-17,2.7166e-17,2.7166e-17,2.7166e-17,2.7166e-17,2.7166e-17,2.7165e-17,2.7165e-17,2.7165e-17,2.7164e-17,2.7163e-17,2.7162e-17,2.7162e-17,2.7161e-17,2.7159e-17,2.7157e-17,2.7156e-17,2.7154e-17,2.7152e-17,2.7149e-17,2.7146e-17,2.7143e-17,2.714e-17,2.7137e-17,2.7133e-17,2.713e-17,2.7126e-17,2.7122e-17,2.7117e-17,2.7113e-17,2.7109e-17,2.7104e-17,2.71e-17,2.7095e-17,2.709e-17,2.7085e-17,2.708e-17,2.7075e-17,2.7069e-17,2.7063e-17,2.7058e-17,2.7052e-17,2.7045e-17,2.7039e-17,2.7032e-17,2.7026e-17,2.7018e-17,2.7011e-17,2.7004e-17,2.6996e-17,2.6988e-17,2.698e-17,2.6971e-17,2.6962e-17,2.6953e-17,2.6944e-17,2.6935e-17,2.6924e-17,2.6913e-17,2.6901e-17,2.6889e-17,2.6876e-17,2.6863e-17,2.6849e-17,2.6835e-17,2.6821e-17,2.6805e-17,2.6789e-17,2.6773e-17,2.6757e-17,2.6741e-17,2.6722e-17,2.6703e-17,2.6684e-17,2.6665e-17,2.6646e-17,2.6625e-17,2.6602e-17,2.658e-17,2.6558e-17,2.6535e-17,2.6511e-17,2.6486e-17,2.646e-17,2.6435e-17,2.641e-17,2.6382e-17,2.6353e-17,2.6323e-17,2.6295e-17,2.6267e-17,2.6237e-17,2.6203e-17,2.6168e-17,2.6135e-17,2.6103e-17,2.607e-17,2.6035e-17,2.5995e-17,2.5956e-17,2.5918e-17,2.5881e-17,2.5842e-17,2.5798e-17,2.5753e-17,2.571e-17,2.5669e-17,2.5626e-17,2.5577e-17,2.5554e-17,2.554e-17,2.5531e-17,2.5525e-17,2.5516e-17,2.5501e-17,2.5484e-17,2.5469e-17,2.5461e-17,2.5455e-17,2.5445e-17,2.543e-17,2.5412e-17,2.5398e-17,2.5389e-17,2.5382e-17,2.5371e-17,2.5354e-17,2.5336e-17,2.5323e-17,2.5315e-17,2.5307e-17,2.5294e-17,2.5275e-17,2.5256e-17,2.5242e-17,2.5236e-17,2.5229e-17,2.5216e-17,2.5195e-17,2.5173e-17,2.5157e-17,2.5149e-17,2.5145e-17,2.5134e-17,2.5113e-17,2.5089e-17,2.507e-17,2.506e-17,2.5055e-17,2.5045e-17,2.5026e-17,2.5002e-17,2.4981e-17,2.497e-17,2.4964e-17,2.4934e-17,2.4893e-17,2.4844e-17,2.4797e-17,2.4761e-17,2.4735e-17,2.4708e-17,2.4668e-17,2.4614e-17,2.4558e-17,2.4511e-17,2.448e-17,2.4455e-17,2.4422e-17,2.437e-17,2.4308e-17,2.4249e-17,2.4206e-17,2.4177e-17,2.4147e-17,2.4102e-17,2.404e-17,2.3972e-17,2.3916e-17,2.3878e-17,2.3848e-17,2.3809e-17,2.375e-17,2.3677e-17,2.3606e-17,2.3554e-17,2.3521e-17,2.3489e-17,2.3441e-17,2.3369e-17,2.3286e-17,2.321e-17,2.3159e-17,2.3126e-17,2.3092e-17,2.3038e-17,2.2979e-17,2.2918e-17,2.2869e-17,2.2846e-17,2.284e-17,2.2833e-17,2.2808e-17,2.2759e-17,2.2699e-17,2.2651e-17,2.2626e-17,2.2617e-17,2.2607e-17,2.2581e-17,2.253e-17,2.2469e-17,2.2419e-17,2.2392e-17,2.2382e-17,2.2375e-17,2.2353e-17,2.2302e-17,2.2236e-17,2.2179e-17,2.2142e-17,2.2127e-17,2.2125e-17,2.2112e-17,2.2071e-17,2.2008e-17,2.1942e-17,2.1889e-17,2.1859e-17,2.1853e-17,2.1849e-17,2.1822e-17,2.177e-17,2.1704e-17,2.1665e-17,2.1649e-17,2.1662e-17,2.1688e-17,2.1704e-17,2.1696e-17,2.1667e-17,2.1628e-17,2.1598e-17,2.1596e-17,2.1619e-17,2.1646e-17,2.1656e-17,2.1639e-17,2.16e-17,2.1555e-17,2.1531e-17,2.1538e-17,2.1568e-17,2.1599e-17,2.1609e-17,2.1588e-17,2.1541e-17,2.1493e-17,2.1468e-17,2.1476e-17,2.1508e-17,2.1542e-17,2.1554e-17,2.1533e-17,2.1487e-17,2.1439e-17,2.1413e-17,2.1418e-17,2.1448e-17,2.1465e-17,2.1462e-17,2.1425e-17,2.1364e-17,2.1301e-17,2.1255e-17,2.1243e-17,2.1259e-17,2.1279e-17,2.1277e-17,2.1244e-17,2.1184e-17,2.1114e-17,2.1058e-17,2.1039e-17,2.1054e-17,2.1078e-17,2.1086e-17,2.1067e-17,2.1017e-17,2.0944e-17,2.0873e-17,2.0835e-17,2.0837e-17,2.0856e-17,2.0872e-17,2.0871e-17,2.0844e-17,2.0783e-17,2.0705e-17,2.0648e-17,2.0632e-17,2.0654e-17,2.0682e-17,2.0707e-17,2.0721e-17,2.0702e-17,2.065e-17,2.0598e-17,2.0582e-17,2.0594e-17,2.0613e-17,2.0636e-17,2.0663e-17,2.0673e-17,2.0642e-17,2.0584e-17,2.0542e-17,2.0536e-17,2.0544e-17,2.0556e-17,2.0582e-17,2.0618e-17,2.063e-17,2.0596e-17,2.054e-17,2.0507e-17,2.0496e-17,2.049e-17,2.0491e-17,2.0514e-17,2.0557e-17,2.0479e-17,2.0345e-17,2.0188e-17,2.005e-17,1.9934e-17,1.9802e-17,1.9658e-17,1.9552e-17,1.9472e-17,1.9382e-17,1.9247e-17,1.9079e-17,1.8937e-17,1.8801e-17,1.865e-17,1.8474e-17,1.8325e-17,1.8226e-17,1.8122e-17,1.7979e-17,1.7793e-17,1.7628e-17,1.7493e-17,1.7331e-17,1.7122e-17,1.6925e-17,1.6777e-17,1.6645e-17,1.6653e-17,1.6644e-17,1.6625e-17,1.664e-17,1.6654e-17,1.6622e-17,1.6567e-17,1.6554e-17,1.6571e-17,1.6592e-17,1.6603e-17,1.6605e-17,1.6612e-17,1.6636e-17,1.6628e-17,1.6587e-17,1.6543e-17,1.6535e-17,1.654e-17,1.654e-17,1.6549e-17,1.6566e-17,1.659e-17,1.6615e-17,1.6612e-17,1.6555e-17,1.6528e-17,1.6521e-17,1.651e-17,1.6495e-17,1.6504e-17,1.6526e-17,1.6563e-17,1.6591e-17,1.6581e-17,1.6539e-17,1.6511e-17,1.6505e-17,1.6488e-17,1.6458e-17,1.6457e-17,1.6476e-17,1.6525e-17,1.6565e-17,1.6567e-17,1.6538e-17,1.6511e-17,1.6503e-17,1.6483e-17,1.6438e-17,1.6407e-17,1.6403e-17,1.6433e-17,1.6473e-17,1.6481e-17,1.6455e-17,1.6417e-17,1.6396e-17,1.6378e-17,1.632e-17,1.6243e-17,1.6201e-17,1.6204e-17,1.624e-17,1.6273e-17,1.6262e-17,1.623e-17,1.6208e-17,1.6205e-17,1.6189e-17,1.6117e-17,1.603e-17,1.5988e-17,1.5989e-17,1.6037e-17,1.6064e-17,1.6041e-17,1.6002e-17,1.5986e-17,1.5994e-17,1.5991e-17,1.5923e-17,1.5821e-17,1.5762e-17,1.5756e-17,1.5802e-17,1.5843e-17,1.5824e-17,1.5773e-17,1.575e-17,1.5758e-17,1.5782e-17,1.5751e-17,1.5647e-17,1.5547e-17,1.5512e-17,1.5537e-17} };

    const static float base_density = 2.3864e+16f;//km^-3

    const static float Hp = 30000e-3f; // Continental value

    const static float pb = 2e6f; // background value in km^-3

    template <typename Float, typename UInt32, typename Mask, typename Spectrum, typename Wavelength>
    class PolarArticAerosol : public GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength> {
        const float pb_bd;

    public:
        PolarArticAerosol() : GlobalAerosolModel<Float, UInt32, Mask, Spectrum, Wavelength>(tabulatedValues),
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
#endif //_AEROSOL_POLAR_ARTIC