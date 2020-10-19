#ifndef _AEROSOL_DESERT_DUST
#define _AEROSOL_DESERT_DUST

#include <vector>
#include <mitsuba/atmosphere/GlobalAtmosphericAerosol.h>
#include <mitsuba/atmosphere/Utils.h>

namespace desertDustAerosol {

	// row[0] = wavelength
	// row[1] = cross section absorption coefficient
	// row[2] = cross section scattering coefficient
    template <typename Float>
	const static Float tabulatedValues[3][1001] =
	{ { 100,100.23,100.46,100.69,100.93,101.16,101.39,101.62,101.86,102.09,102.33,102.57,102.8,103.04,103.28,103.51,103.75,103.99,104.23,104.47,104.71,104.95,105.2,105.44,105.68,105.93,106.17,106.41,106.66,106.91,107.15,107.4,107.65,107.89,108.14,108.39,108.64,108.89,109.14,109.4,109.65,109.9,110.15,110.41,110.66,110.92,111.17,111.43,111.69,111.94,112.2,112.46,112.72,112.98,113.24,113.5,113.76,114.03,114.29,114.55,114.82,115.08,115.35,115.61,115.88,116.14,116.41,116.68,116.95,117.22,117.49,117.76,118.03,118.3,118.58,118.85,119.12,119.4,119.67,119.95,120.23,120.5,120.78,121.06,121.34,121.62,121.9,122.18,122.46,122.74,123.03,123.31,123.59,123.88,124.17,124.45,124.74,125.03,125.31,125.6,125.89,126.18,126.47,126.77,127.06,127.35,127.64,127.94,128.23,128.53,128.82,129.12,129.42,129.72,130.02,130.32,130.62,130.92,131.22,131.52,131.83,132.13,132.43,132.74,133.05,133.35,133.66,133.97,134.28,134.59,134.9,135.21,135.52,135.83,136.14,136.46,136.77,137.09,137.4,137.72,138.04,138.36,138.68,139,139.32,139.64,139.96,140.28,140.6,140.93,141.25,141.58,141.91,142.23,142.56,142.89,143.22,143.55,143.88,144.21,144.54,144.88,145.21,145.55,145.88,146.22,146.55,146.89,147.23,147.57,147.91,148.25,148.59,148.94,149.28,149.62,149.97,150.31,150.66,151.01,151.36,151.71,152.05,152.41,152.76,153.11,153.46,153.82,154.17,154.53,154.88,155.24,155.6,155.96,156.31,156.68,157.04,157.4,157.76,158.12,158.49,158.85,159.22,159.59,159.96,160.32,160.69,161.06,161.44,161.81,162.18,162.55,162.93,163.31,163.68,164.06,164.44,164.82,165.2,165.58,165.96,166.34,166.72,167.11,167.49,167.88,168.27,168.66,169.04,169.43,169.82,170.22,170.61,171,171.4,171.79,172.19,172.58,172.98,173.38,173.78,174.18,174.58,174.98,175.39,175.79,176.2,176.6,177.01,177.42,177.83,178.24,178.65,179.06,179.47,179.89,180.3,180.72,181.13,181.55,181.97,182.39,182.81,183.23,183.65,184.08,184.5,184.93,185.35,185.78,186.21,186.64,187.07,187.5,187.93,188.36,188.8,189.23,189.67,190.11,190.55,190.99,191.43,191.87,192.31,192.75,193.2,193.64,194.09,194.54,194.98,195.43,195.88,196.34,196.79,197.24,197.7,198.15,198.61,199.07,199.53,199.99,200.45,200.91,201.37,201.84,202.3,202.77,203.24,203.7,204.17,204.64,205.12,205.59,206.06,206.54,207.01,207.49,207.97,208.45,208.93,209.41,209.89,210.38,210.86,211.35,211.84,212.32,212.81,213.3,213.8,214.29,214.78,215.28,215.77,216.27,216.77,217.27,217.77,218.27,218.78,219.28,219.79,220.29,220.8,221.31,221.82,222.33,222.84,223.36,223.87,224.39,224.91,225.42,225.94,226.46,226.99,227.51,228.03,228.56,229.09,229.61,230.14,230.67,231.21,231.74,232.27,232.81,233.35,233.88,234.42,234.96,235.5,236.05,236.59,237.14,237.68,238.23,238.78,239.33,239.88,240.44,240.99,241.55,242.1,242.66,243.22,243.78,244.34,244.91,245.47,246.04,246.6,247.17,247.74,248.31,248.89,249.46,250.03,250.61,251.19,251.77,252.35,252.93,253.51,254.1,254.68,255.27,255.86,256.45,257.04,257.63,258.23,258.82,259.42,260.02,260.62,261.22,261.82,262.42,263.03,263.63,264.24,264.85,265.46,266.07,266.69,267.3,267.92,268.53,269.15,269.77,270.4,271.02,271.64,272.27,272.9,273.53,274.16,274.79,275.42,276.06,276.69,277.33,277.97,278.61,279.25,279.9,280.54,281.19,281.84,282.49,283.14,283.79,284.45,285.1,285.76,286.42,287.08,287.74,288.4,289.07,289.73,290.4,291.07,291.74,292.42,293.09,293.76,294.44,295.12,295.8,296.48,297.17,297.85,298.54,299.23,299.92,300.61,301.3,302,302.69,303.39,304.09,304.79,305.49,306.2,306.9,307.61,308.32,309.03,309.74,310.46,311.17,311.89,312.61,313.33,314.05,314.77,315.5,316.23,316.96,317.69,318.42,319.15,319.89,320.63,321.37,322.11,322.85,323.59,324.34,325.09,325.84,326.59,327.34,328.1,328.85,329.61,330.37,331.13,331.89,332.66,333.43,334.19,334.97,335.74,336.51,337.29,338.06,338.84,339.63,340.41,341.19,341.98,342.77,343.56,344.35,345.14,345.94,346.74,347.54,348.34,349.14,349.95,350.75,351.56,352.37,353.18,354,354.81,355.63,356.45,357.27,358.1,358.92,359.75,360.58,361.41,362.24,363.08,363.92,364.75,365.59,366.44,367.28,368.13,368.98,369.83,370.68,371.54,372.39,373.25,374.11,374.97,375.84,376.7,377.57,378.44,379.32,380.19,381.07,381.94,382.82,383.71,384.59,385.48,386.37,387.26,388.15,389.05,389.94,390.84,391.74,392.64,393.55,394.46,395.37,396.28,397.19,398.11,399.02,399.94,400.87,401.79,402.72,403.65,404.58,405.51,406.44,407.38,408.32,409.26,410.2,411.15,412.1,413.05,414,414.95,415.91,416.87,417.83,418.79,419.76,420.73,421.7,422.67,423.64,424.62,425.6,426.58,427.56,428.55,429.54,430.53,431.52,432.51,433.51,434.51,435.51,436.52,437.52,438.53,439.54,440.55,441.57,442.59,443.61,444.63,445.66,446.68,447.71,448.75,449.78,450.82,451.86,452.9,453.94,454.99,456.04,457.09,458.14,459.2,460.26,461.32,462.38,463.45,464.52,465.59,466.66,467.74,468.81,469.89,470.98,472.06,473.15,474.24,475.34,476.43,477.53,478.63,479.73,480.84,481.95,483.06,484.17,485.29,486.41,487.53,488.65,489.78,490.91,492.04,493.17,494.31,495.45,496.59,497.74,498.88,500.03,501.19,502.34,503.5,504.66,505.82,506.99,508.16,509.33,510.51,511.68,512.86,514.04,515.23,516.42,517.61,518.8,520,521.19,522.4,523.6,524.81,526.02,527.23,528.45,529.66,530.88,532.11,533.33,534.56,535.8,537.03,538.27,539.51,540.75,542,543.25,544.5,545.76,547.02,548.28,549.54,550.81,552.08,553.35,554.63,555.9,557.19,558.47,559.76,561.05,562.34,563.64,564.94,566.24,567.54,568.85,570.16,571.48,572.8,574.12,575.44,576.77,578.1,579.43,580.76,582.1,583.45,584.79,586.14,587.49,588.84,590.2,591.56,592.93,594.29,595.66,597.04,598.41,599.79,601.17,602.56,603.95,605.34,606.74,608.13,609.54,610.94,612.35,613.76,615.18,616.6,618.02,619.44,620.87,622.3,623.73,625.17,626.61,628.06,629.51,630.96,632.41,633.87,635.33,636.8,638.26,639.73,641.21,642.69,644.17,645.65,647.14,648.63,650.13,651.63,653.13,654.64,656.15,657.66,659.17,660.69,662.22,663.74,665.27,666.81,668.34,669.88,671.43,672.98,674.53,676.08,677.64,679.2,680.77,682.34,683.91,685.49,687.07,688.65,690.24,691.83,693.43,695.02,696.63,698.23,699.84,701.46,703.07,704.69,706.32,707.95,709.58,711.21,712.85,714.5,716.14,717.79,719.45,721.11,722.77,724.44,726.11,727.78,729.46,731.14,732.82,734.51,736.21,737.9,739.61,741.31,743.02,744.73,746.45,748.17,749.89,751.62,753.36,755.09,756.83,758.58,760.33,762.08,763.84,765.6,767.36,769.13,770.9,772.68,774.46,776.25,778.04,779.83,781.63,783.43,785.24,787.05,788.86,790.68,792.5,794.33,796.16,797.99,799.83,801.68,803.53,805.38,807.24,809.1,810.96,812.83,814.7,816.58,818.46,820.35,822.24,824.14,826.04,827.94,829.85,831.76,833.68,835.6,837.53,839.46,841.4,843.33,845.28,847.23,849.18,851.14,853.1,855.07,857.04,859.01,860.99,862.98,864.97,866.96,868.96,870.96,872.97,874.98,877,879.02,881.05,883.08,885.12,887.16,889.2,891.25,893.31,895.36,897.43,899.5,901.57,903.65,905.73,907.82,909.91,912.01,914.11,916.22,918.33,920.45,922.57,924.7,926.83,928.97,931.11,933.25,935.41,937.56,939.72,941.89,944.06,946.24,948.42,950.6,952.8,954.99,957.19,959.4,961.61,963.83,966.05,968.28,970.51,972.75,974.99,977.24,979.49,981.75,984.01,986.28,988.55,990.83,993.12,995.41,997.7,1000 },
	{ 4.0763e-16,4.0763e-16,4.0764e-16,4.0765e-16,4.0765e-16,4.0766e-16,4.0767e-16,4.0767e-16,4.0768e-16,4.0769e-16,4.077e-16,4.077e-16,4.0771e-16,4.0772e-16,4.0772e-16,4.0773e-16,4.0774e-16,4.0775e-16,4.0775e-16,4.0776e-16,4.0777e-16,4.0777e-16,4.0778e-16,4.0779e-16,4.078e-16,4.078e-16,4.0781e-16,4.0782e-16,4.0783e-16,4.0783e-16,4.0784e-16,4.0785e-16,4.0785e-16,4.0786e-16,4.0787e-16,4.0788e-16,4.0788e-16,4.0789e-16,4.079e-16,4.0791e-16,4.0791e-16,4.0792e-16,4.0793e-16,4.0793e-16,4.0794e-16,4.0795e-16,4.0796e-16,4.0796e-16,4.0797e-16,4.0798e-16,4.0799e-16,4.0799e-16,4.08e-16,4.0801e-16,4.0802e-16,4.0802e-16,4.0803e-16,4.0804e-16,4.0805e-16,4.0805e-16,4.0806e-16,4.0807e-16,4.0808e-16,4.0808e-16,4.0809e-16,4.081e-16,4.0811e-16,4.0811e-16,4.0812e-16,4.0813e-16,4.0814e-16,4.0814e-16,4.0815e-16,4.0816e-16,4.0817e-16,4.0817e-16,4.0818e-16,4.0819e-16,4.082e-16,4.082e-16,4.0821e-16,4.0822e-16,4.0823e-16,4.0823e-16,4.0824e-16,4.0825e-16,4.0826e-16,4.0826e-16,4.0827e-16,4.0828e-16,4.0829e-16,4.083e-16,4.083e-16,4.0831e-16,4.0832e-16,4.0833e-16,4.0833e-16,4.0834e-16,4.0835e-16,4.0836e-16,4.0836e-16,4.0837e-16,4.0838e-16,4.0839e-16,4.084e-16,4.084e-16,4.0841e-16,4.0842e-16,4.0843e-16,4.0843e-16,4.0844e-16,4.0845e-16,4.0846e-16,4.0847e-16,4.0847e-16,4.0848e-16,4.0849e-16,4.085e-16,4.0851e-16,4.0851e-16,4.0852e-16,4.0853e-16,4.0854e-16,4.0855e-16,4.0855e-16,4.0856e-16,4.0857e-16,4.0858e-16,4.0858e-16,4.0859e-16,4.086e-16,4.0861e-16,4.0862e-16,4.0862e-16,4.0863e-16,4.0864e-16,4.0865e-16,4.0866e-16,4.0866e-16,4.0867e-16,4.0868e-16,4.0869e-16,4.087e-16,4.0871e-16,4.0871e-16,4.0872e-16,4.0873e-16,4.0874e-16,4.0875e-16,4.0875e-16,4.0876e-16,4.0877e-16,4.0878e-16,4.0879e-16,4.0879e-16,4.088e-16,4.0881e-16,4.0882e-16,4.0883e-16,4.0884e-16,4.0884e-16,4.0885e-16,4.0886e-16,4.0887e-16,4.0888e-16,4.0888e-16,4.0889e-16,4.089e-16,4.0891e-16,4.0892e-16,4.0893e-16,4.0893e-16,4.0894e-16,4.0895e-16,4.0896e-16,4.0897e-16,4.0898e-16,4.0898e-16,4.0899e-16,4.09e-16,4.0901e-16,4.0902e-16,4.0903e-16,4.0903e-16,4.0904e-16,4.0905e-16,4.0906e-16,4.0907e-16,4.0908e-16,4.0908e-16,4.0909e-16,4.091e-16,4.0911e-16,4.0912e-16,4.0913e-16,4.0914e-16,4.0914e-16,4.0915e-16,4.0916e-16,4.0917e-16,4.0918e-16,4.0919e-16,4.0919e-16,4.092e-16,4.0921e-16,4.0922e-16,4.0923e-16,4.0924e-16,4.0925e-16,4.0925e-16,4.0926e-16,4.0927e-16,4.0928e-16,4.0929e-16,4.093e-16,4.0931e-16,4.0931e-16,4.0932e-16,4.0933e-16,4.0934e-16,4.0935e-16,4.0936e-16,4.0937e-16,4.0937e-16,4.0938e-16,4.0939e-16,4.094e-16,4.0941e-16,4.0942e-16,4.0943e-16,4.0944e-16,4.0944e-16,4.0945e-16,4.0946e-16,4.0947e-16,4.0948e-16,4.0949e-16,4.095e-16,4.0951e-16,4.0951e-16,4.0952e-16,4.0953e-16,4.0954e-16,4.0955e-16,4.0956e-16,4.0957e-16,4.0958e-16,4.0959e-16,4.0959e-16,4.096e-16,4.0961e-16,4.0962e-16,4.0963e-16,4.0964e-16,4.0965e-16,4.0966e-16,4.0967e-16,4.0967e-16,4.0968e-16,4.0969e-16,4.097e-16,4.0971e-16,4.0972e-16,4.0973e-16,4.0974e-16,4.0975e-16,4.0975e-16,4.0976e-16,4.0977e-16,4.0978e-16,4.0979e-16,4.098e-16,4.0981e-16,4.0982e-16,4.0983e-16,4.0984e-16,4.0985e-16,4.0985e-16,4.0986e-16,4.0987e-16,4.0988e-16,4.0989e-16,4.099e-16,4.0991e-16,4.0992e-16,4.0993e-16,4.0994e-16,4.0995e-16,4.0995e-16,4.0996e-16,4.0997e-16,4.0998e-16,4.0999e-16,4.1e-16,4.1001e-16,4.1002e-16,4.1003e-16,4.1004e-16,4.1005e-16,4.1006e-16,4.1007e-16,4.1007e-16,4.1008e-16,4.1009e-16,4.101e-16,4.1011e-16,4.1012e-16,4.1013e-16,4.1014e-16,4.1015e-16,4.1016e-16,4.1017e-16,4.1018e-16,4.1019e-16,4.102e-16,4.1021e-16,4.1021e-16,4.1022e-16,4.1023e-16,4.1024e-16,4.1025e-16,4.1026e-16,4.1027e-16,4.1028e-16,4.1029e-16,4.103e-16,4.1031e-16,4.1032e-16,4.1033e-16,4.1034e-16,4.1035e-16,4.1036e-16,4.1037e-16,4.1038e-16,4.1038e-16,4.1039e-16,4.104e-16,4.1041e-16,4.1042e-16,4.1043e-16,4.1044e-16,4.1045e-16,4.1046e-16,4.1047e-16,4.1048e-16,4.1049e-16,4.105e-16,4.1051e-16,4.1052e-16,4.1053e-16,4.1054e-16,4.1055e-16,4.1056e-16,4.1057e-16,4.1058e-16,4.1059e-16,4.106e-16,4.1061e-16,4.1062e-16,4.1063e-16,4.1064e-16,4.1064e-16,4.1065e-16,4.1066e-16,4.1067e-16,4.1068e-16,4.1069e-16,4.107e-16,4.1071e-16,4.1072e-16,4.1073e-16,4.1074e-16,4.1075e-16,4.1076e-16,4.1077e-16,4.1078e-16,4.1079e-16,4.108e-16,4.1081e-16,4.1082e-16,4.1083e-16,4.1084e-16,4.1085e-16,4.1086e-16,4.1087e-16,4.1088e-16,4.1089e-16,4.109e-16,4.1091e-16,4.1092e-16,4.1093e-16,4.1094e-16,4.1095e-16,4.1096e-16,4.1097e-16,4.1098e-16,4.1099e-16,4.11e-16,4.1101e-16,4.1102e-16,4.1103e-16,4.1104e-16,4.1105e-16,4.1106e-16,4.1107e-16,4.1108e-16,4.1109e-16,4.111e-16,4.1111e-16,4.1112e-16,4.1113e-16,4.1114e-16,4.1115e-16,4.1116e-16,4.1117e-16,4.1118e-16,4.1119e-16,4.112e-16,4.1121e-16,4.1122e-16,4.1123e-16,4.1124e-16,4.1125e-16,4.1126e-16,4.1127e-16,4.1128e-16,4.1129e-16,4.1131e-16,4.1132e-16,4.1133e-16,4.1134e-16,4.1135e-16,4.1136e-16,4.1137e-16,4.1138e-16,4.1139e-16,4.114e-16,4.1141e-16,4.1142e-16,4.1143e-16,4.1144e-16,4.1145e-16,4.1146e-16,4.1147e-16,4.1148e-16,4.1149e-16,4.115e-16,4.1151e-16,4.1152e-16,4.1153e-16,4.1154e-16,4.1155e-16,4.1156e-16,4.1157e-16,4.1159e-16,4.116e-16,4.1161e-16,4.1162e-16,4.1163e-16,4.1164e-16,4.1165e-16,4.1166e-16,4.1167e-16,4.1168e-16,4.1169e-16,4.117e-16,4.1171e-16,4.1172e-16,4.1173e-16,4.1174e-16,4.1175e-16,4.1177e-16,4.1178e-16,4.1179e-16,4.118e-16,4.1181e-16,4.1182e-16,4.1183e-16,4.1184e-16,4.1185e-16,4.1186e-16,4.1187e-16,4.1188e-16,4.1189e-16,4.119e-16,4.1191e-16,4.1193e-16,4.1194e-16,4.1195e-16,4.1196e-16,4.1197e-16,4.1198e-16,4.1199e-16,4.12e-16,4.1201e-16,4.1202e-16,4.1203e-16,4.1204e-16,4.1205e-16,4.1206e-16,4.1208e-16,4.1209e-16,4.121e-16,4.1211e-16,4.1212e-16,4.1213e-16,4.1214e-16,4.1215e-16,4.1216e-16,4.1217e-16,4.1218e-16,4.122e-16,4.1221e-16,4.1222e-16,4.1223e-16,4.1224e-16,4.1225e-16,4.1227e-16,4.1228e-16,4.1229e-16,4.123e-16,4.1231e-16,4.1232e-16,4.1233e-16,4.1234e-16,4.1235e-16,4.1236e-16,4.1237e-16,4.1239e-16,4.124e-16,4.1241e-16,4.1243e-16,4.1244e-16,4.1244e-16,4.1245e-16,4.1246e-16,4.1247e-16,4.1248e-16,4.1249e-16,4.1251e-16,4.1252e-16,4.1254e-16,4.1255e-16,4.1257e-16,4.1258e-16,4.1259e-16,4.1259e-16,4.126e-16,4.1261e-16,4.1262e-16,4.1262e-16,4.1264e-16,4.1265e-16,4.1267e-16,4.1269e-16,4.1271e-16,4.1272e-16,4.1274e-16,4.1275e-16,4.1275e-16,4.1276e-16,4.1276e-16,4.1277e-16,4.1277e-16,4.1278e-16,4.128e-16,4.1282e-16,4.1284e-16,4.1286e-16,4.1289e-16,4.1291e-16,4.1292e-16,4.1293e-16,4.1294e-16,4.1294e-16,4.1293e-16,4.1293e-16,4.1294e-16,4.1294e-16,4.1296e-16,4.1299e-16,4.1302e-16,4.1305e-16,4.1309e-16,4.1312e-16,4.1314e-16,4.1315e-16,4.1316e-16,4.1315e-16,4.1314e-16,4.1313e-16,4.1312e-16,4.1312e-16,4.1314e-16,4.1317e-16,4.1321e-16,4.1325e-16,4.1331e-16,4.1336e-16,4.134e-16,4.1344e-16,4.1345e-16,4.1345e-16,4.1343e-16,4.1341e-16,4.1339e-16,4.1337e-16,4.1337e-16,4.1338e-16,4.1342e-16,4.1349e-16,4.1357e-16,4.1366e-16,4.1375e-16,4.1384e-16,4.139e-16,4.1394e-16,4.1395e-16,4.1393e-16,4.1389e-16,4.1385e-16,4.138e-16,4.1378e-16,4.1379e-16,4.1383e-16,4.1392e-16,4.1405e-16,4.1421e-16,4.1438e-16,4.1455e-16,4.1471e-16,4.1483e-16,4.149e-16,4.1493e-16,4.1491e-16,4.1485e-16,4.1478e-16,4.1472e-16,4.1468e-16,4.147e-16,4.1478e-16,4.1494e-16,4.1518e-16,4.1547e-16,4.1581e-16,4.1616e-16,4.1649e-16,4.1678e-16,4.17e-16,4.1713e-16,4.1719e-16,4.1717e-16,4.1711e-16,4.1702e-16,4.1686e-16,4.1674e-16,4.167e-16,4.1678e-16,4.1697e-16,4.1725e-16,4.1762e-16,4.1802e-16,4.1842e-16,4.1877e-16,4.1905e-16,4.1922e-16,4.1927e-16,4.192e-16,4.1904e-16,4.1881e-16,4.1854e-16,4.183e-16,4.1811e-16,4.1801e-16,4.1803e-16,4.1819e-16,4.1847e-16,4.1886e-16,4.1932e-16,4.1982e-16,4.203e-16,4.2072e-16,4.2104e-16,4.2124e-16,4.2129e-16,4.2122e-16,4.2102e-16,4.2075e-16,4.2044e-16,4.2013e-16,4.1989e-16,4.1975e-16,4.1975e-16,4.1989e-16,4.2019e-16,4.2062e-16,4.2115e-16,4.2175e-16,4.2236e-16,4.2292e-16,4.2355e-16,4.2405e-16,4.2441e-16,4.246e-16,4.2464e-16,4.2455e-16,4.2437e-16,4.2414e-16,4.2393e-16,4.2379e-16,4.2379e-16,4.2396e-16,4.2434e-16,4.2493e-16,4.2572e-16,4.2668e-16,4.2776e-16,4.289e-16,4.3002e-16,4.3108e-16,4.3199e-16,4.3272e-16,4.3323e-16,4.3352e-16,4.336e-16,4.3352e-16,4.3333e-16,4.331e-16,4.3292e-16,4.3287e-16,4.3301e-16,4.3341e-16,4.3409e-16,4.3508e-16,4.3636e-16,4.3788e-16,4.3959e-16,4.414e-16,4.4323e-16,4.4498e-16,4.4656e-16,4.4763e-16,4.4823e-16,4.4852e-16,4.485e-16,4.4821e-16,4.4773e-16,4.4714e-16,4.4654e-16,4.4603e-16,4.457e-16,4.4564e-16,4.4589e-16,4.4651e-16,4.4748e-16,4.4879e-16,4.5038e-16,4.5217e-16,4.5406e-16,4.5596e-16,4.5775e-16,4.5934e-16,4.6065e-16,4.6162e-16,4.6222e-16,4.6244e-16,4.6233e-16,4.6193e-16,4.6134e-16,4.6066e-16,4.6e-16,4.5946e-16,4.5916e-16,4.5917e-16,4.5957e-16,4.6039e-16,4.6163e-16,4.6327e-16,4.6526e-16,4.6699e-16,4.6877e-16,4.7056e-16,4.7224e-16,4.737e-16,4.7486e-16,4.7563e-16,4.7599e-16,4.7593e-16,4.7545e-16,4.7463e-16,4.7352e-16,4.7222e-16,4.7086e-16,4.6954e-16,4.6836e-16,4.6744e-16,4.6686e-16,4.6667e-16,4.6691e-16,4.6758e-16,4.6865e-16,4.7007e-16,4.7175e-16,4.7361e-16,4.7553e-16,4.7739e-16,4.7911e-16,4.8056e-16,4.8167e-16,4.8239e-16,4.8269e-16,4.8254e-16,4.82e-16,4.8113e-16,4.8029e-16,4.7926e-16,4.7814e-16,4.7703e-16,4.7604e-16,4.7529e-16,4.7486e-16,4.7481e-16,4.752e-16,4.7604e-16,4.7732e-16,4.7902e-16,4.8106e-16,4.8337e-16,4.8585e-16,4.8838e-16,4.9087e-16,4.9319e-16,4.9525e-16,4.9696e-16,4.9826e-16,4.9911e-16,4.995e-16,4.9944e-16,4.9897e-16,4.9815e-16,4.9707e-16,4.9585e-16,4.9456e-16,4.9335e-16,4.9231e-16,4.9155e-16,4.9076e-16,4.9036e-16,4.9043e-16,4.9098e-16,4.9201e-16,4.9347e-16,4.9529e-16,4.974e-16,4.997e-16,5.0209e-16,5.0444e-16,5.0664e-16,5.0861e-16,5.1025e-16,5.1149e-16,5.1226e-16,5.1257e-16,5.1239e-16,5.1176e-16,5.1071e-16,5.0932e-16,5.0768e-16,5.0586e-16,5.04e-16,5.0219e-16,5.0053e-16,4.9913e-16,4.9806e-16,4.9741e-16,4.9721e-16,5.0008e-16,5.0367e-16,5.079e-16,5.127e-16,5.1805e-16,5.2389e-16,5.3015e-16,5.367e-16,5.4345e-16,5.5029e-16,5.5709e-16,5.6373e-16,5.7012e-16,5.762e-16,5.8186e-16,5.8707e-16,5.9182e-16,5.9617e-16,6.0011e-16,6.0371e-16,6.0709e-16,6.104e-16,6.1374e-16,6.1724e-16,6.2113e-16,6.2554e-16,6.3062e-16,6.3645e-16,6.3616e-16,6.3602e-16,6.3669e-16,6.3804e-16,6.4011e-16,6.4279e-16,6.4607e-16,6.4971e-16,6.5357e-16,6.5751e-16,6.6151e-16,6.6531e-16,6.6872e-16,6.716e-16,6.7401e-16,6.758e-16,6.7681e-16,6.7699e-16,6.7653e-16,6.754e-16,6.7353e-16,6.7102e-16,6.6813e-16,6.6495e-16,6.6146e-16,6.5782e-16,6.5435e-16,6.5117e-16,6.4822e-16,6.4567e-16,6.4377e-16,6.426e-16,6.4201e-16,6.4208e-16,6.4298e-16,6.4468e-16,6.4691e-16,6.4965e-16,6.5295e-16,6.5675e-16,6.6073e-16,6.6475e-16,6.6875e-16,6.7272e-16,6.7643e-16,6.7967e-16,6.8238e-16,6.8454e-16,6.8611e-16,6.8695e-16,6.8709e-16,6.8671e-16,6.8577e-16,6.8417e-16,6.8208e-16,6.7945e-16,6.7646e-16,6.7321e-16,6.6988e-16,6.6656e-16,6.633e-16,6.6031e-16,6.5764e-16,6.556e-16,6.5401e-16,6.5296e-16,6.5255e-16,6.53e-16,6.5415e-16,6.5585e-16,6.5822e-16,6.6115e-16,6.6477e-16,6.6872e-16,6.7286e-16,6.7721e-16,6.819e-16,6.865e-16,6.9074e-16,6.9476e-16,6.9841e-16,7.0189e-16,7.0471e-16,7.0665e-16,7.0794e-16,7.0882e-16,7.0912e-16,7.0855e-16,7.0723e-16,7.052e-16,7.0304e-16,7.0045e-16,6.9724e-16,6.9377e-16,6.9013e-16,6.8674e-16,6.8343e-16},
	{ 3.3464e-16,3.3465e-16,3.3466e-16,3.3466e-16,3.3467e-16,3.3468e-16,3.3469e-16,3.3469e-16,3.347e-16,3.3471e-16,3.3472e-16,3.3472e-16,3.3473e-16,3.3474e-16,3.3475e-16,3.3475e-16,3.3476e-16,3.3477e-16,3.3478e-16,3.3478e-16,3.3479e-16,3.348e-16,3.3481e-16,3.3481e-16,3.3482e-16,3.3483e-16,3.3484e-16,3.3484e-16,3.3485e-16,3.3486e-16,3.3487e-16,3.3488e-16,3.3488e-16,3.3489e-16,3.349e-16,3.3491e-16,3.3491e-16,3.3492e-16,3.3493e-16,3.3494e-16,3.3495e-16,3.3495e-16,3.3496e-16,3.3497e-16,3.3498e-16,3.3499e-16,3.3499e-16,3.35e-16,3.3501e-16,3.3502e-16,3.3503e-16,3.3503e-16,3.3504e-16,3.3505e-16,3.3506e-16,3.3507e-16,3.3508e-16,3.3508e-16,3.3509e-16,3.351e-16,3.3511e-16,3.3512e-16,3.3512e-16,3.3513e-16,3.3514e-16,3.3515e-16,3.3516e-16,3.3517e-16,3.3517e-16,3.3518e-16,3.3519e-16,3.352e-16,3.3521e-16,3.3522e-16,3.3523e-16,3.3523e-16,3.3524e-16,3.3525e-16,3.3526e-16,3.3527e-16,3.3528e-16,3.3529e-16,3.3529e-16,3.353e-16,3.3531e-16,3.3532e-16,3.3533e-16,3.3534e-16,3.3535e-16,3.3535e-16,3.3536e-16,3.3537e-16,3.3538e-16,3.3539e-16,3.354e-16,3.3541e-16,3.3542e-16,3.3543e-16,3.3543e-16,3.3544e-16,3.3545e-16,3.3546e-16,3.3547e-16,3.3548e-16,3.3549e-16,3.355e-16,3.3551e-16,3.3552e-16,3.3552e-16,3.3553e-16,3.3554e-16,3.3555e-16,3.3556e-16,3.3557e-16,3.3558e-16,3.3559e-16,3.356e-16,3.3561e-16,3.3562e-16,3.3563e-16,3.3564e-16,3.3564e-16,3.3565e-16,3.3566e-16,3.3567e-16,3.3568e-16,3.3569e-16,3.357e-16,3.3571e-16,3.3572e-16,3.3573e-16,3.3574e-16,3.3575e-16,3.3576e-16,3.3577e-16,3.3578e-16,3.3579e-16,3.358e-16,3.3581e-16,3.3582e-16,3.3583e-16,3.3584e-16,3.3585e-16,3.3585e-16,3.3586e-16,3.3587e-16,3.3588e-16,3.3589e-16,3.359e-16,3.3591e-16,3.3592e-16,3.3593e-16,3.3594e-16,3.3595e-16,3.3596e-16,3.3597e-16,3.3598e-16,3.3599e-16,3.36e-16,3.3601e-16,3.3602e-16,3.3603e-16,3.3604e-16,3.3605e-16,3.3606e-16,3.3607e-16,3.3608e-16,3.361e-16,3.3611e-16,3.3612e-16,3.3613e-16,3.3614e-16,3.3615e-16,3.3616e-16,3.3617e-16,3.3618e-16,3.3619e-16,3.362e-16,3.3621e-16,3.3622e-16,3.3623e-16,3.3624e-16,3.3625e-16,3.3626e-16,3.3627e-16,3.3628e-16,3.3629e-16,3.363e-16,3.3632e-16,3.3633e-16,3.3634e-16,3.3635e-16,3.3636e-16,3.3637e-16,3.3638e-16,3.3639e-16,3.364e-16,3.3641e-16,3.3642e-16,3.3643e-16,3.3645e-16,3.3646e-16,3.3647e-16,3.3648e-16,3.3649e-16,3.365e-16,3.3651e-16,3.3652e-16,3.3653e-16,3.3655e-16,3.3656e-16,3.3657e-16,3.3658e-16,3.3659e-16,3.366e-16,3.3661e-16,3.3662e-16,3.3664e-16,3.3665e-16,3.3666e-16,3.3667e-16,3.3668e-16,3.3669e-16,3.367e-16,3.3672e-16,3.3673e-16,3.3674e-16,3.3675e-16,3.3676e-16,3.3677e-16,3.3679e-16,3.368e-16,3.3681e-16,3.3682e-16,3.3683e-16,3.3684e-16,3.3686e-16,3.3687e-16,3.3688e-16,3.3689e-16,3.369e-16,3.3691e-16,3.3693e-16,3.3694e-16,3.3695e-16,3.3696e-16,3.3697e-16,3.3699e-16,3.37e-16,3.3701e-16,3.3702e-16,3.3704e-16,3.3705e-16,3.3706e-16,3.3707e-16,3.3708e-16,3.371e-16,3.3711e-16,3.3712e-16,3.3713e-16,3.3715e-16,3.3716e-16,3.3717e-16,3.3718e-16,3.372e-16,3.3721e-16,3.3722e-16,3.3723e-16,3.3725e-16,3.3726e-16,3.3727e-16,3.3728e-16,3.373e-16,3.3731e-16,3.3732e-16,3.3733e-16,3.3735e-16,3.3736e-16,3.3737e-16,3.3739e-16,3.374e-16,3.3741e-16,3.3742e-16,3.3744e-16,3.3745e-16,3.3746e-16,3.3748e-16,3.3749e-16,3.375e-16,3.3752e-16,3.3753e-16,3.3754e-16,3.3755e-16,3.3757e-16,3.3758e-16,3.3759e-16,3.3761e-16,3.3762e-16,3.3763e-16,3.3765e-16,3.3766e-16,3.3767e-16,3.3769e-16,3.377e-16,3.3771e-16,3.3773e-16,3.3774e-16,3.3776e-16,3.3777e-16,3.3778e-16,3.378e-16,3.3781e-16,3.3782e-16,3.3784e-16,3.3785e-16,3.3787e-16,3.3788e-16,3.3789e-16,3.3791e-16,3.3792e-16,3.3794e-16,3.3795e-16,3.3796e-16,3.3798e-16,3.3799e-16,3.3801e-16,3.3802e-16,3.3803e-16,3.3805e-16,3.3806e-16,3.3808e-16,3.3809e-16,3.3811e-16,3.3812e-16,3.3813e-16,3.3815e-16,3.3816e-16,3.3818e-16,3.3819e-16,3.3821e-16,3.3822e-16,3.3824e-16,3.3825e-16,3.3826e-16,3.3828e-16,3.3829e-16,3.3831e-16,3.3832e-16,3.3834e-16,3.3835e-16,3.3837e-16,3.3838e-16,3.384e-16,3.3841e-16,3.3843e-16,3.3844e-16,3.3846e-16,3.3847e-16,3.3849e-16,3.385e-16,3.3852e-16,3.3853e-16,3.3855e-16,3.3856e-16,3.3858e-16,3.3859e-16,3.3861e-16,3.3862e-16,3.3864e-16,3.3866e-16,3.3867e-16,3.3869e-16,3.387e-16,3.3872e-16,3.3873e-16,3.3875e-16,3.3876e-16,3.3878e-16,3.388e-16,3.3881e-16,3.3883e-16,3.3884e-16,3.3886e-16,3.3887e-16,3.3889e-16,3.3891e-16,3.3892e-16,3.3894e-16,3.3895e-16,3.3897e-16,3.3899e-16,3.39e-16,3.3902e-16,3.3903e-16,3.3905e-16,3.3907e-16,3.3908e-16,3.391e-16,3.3912e-16,3.3913e-16,3.3915e-16,3.3916e-16,3.3918e-16,3.392e-16,3.3921e-16,3.3923e-16,3.3925e-16,3.3926e-16,3.3928e-16,3.393e-16,3.3931e-16,3.3933e-16,3.3935e-16,3.3936e-16,3.3938e-16,3.394e-16,3.3941e-16,3.3943e-16,3.3945e-16,3.3947e-16,3.3948e-16,3.395e-16,3.3952e-16,3.3953e-16,3.3955e-16,3.3957e-16,3.3959e-16,3.396e-16,3.3962e-16,3.3964e-16,3.3966e-16,3.3967e-16,3.3969e-16,3.3971e-16,3.3973e-16,3.3974e-16,3.3976e-16,3.3978e-16,3.398e-16,3.3981e-16,3.3983e-16,3.3985e-16,3.3987e-16,3.3988e-16,3.399e-16,3.3992e-16,3.3994e-16,3.3996e-16,3.3997e-16,3.3999e-16,3.4001e-16,3.4003e-16,3.4005e-16,3.4007e-16,3.4008e-16,3.401e-16,3.4012e-16,3.4014e-16,3.4016e-16,3.4018e-16,3.4019e-16,3.4021e-16,3.4023e-16,3.4025e-16,3.4027e-16,3.4029e-16,3.4031e-16,3.4032e-16,3.4034e-16,3.4036e-16,3.4038e-16,3.404e-16,3.4042e-16,3.4044e-16,3.4046e-16,3.4047e-16,3.4049e-16,3.4051e-16,3.4053e-16,3.4055e-16,3.4057e-16,3.4059e-16,3.4061e-16,3.4063e-16,3.4065e-16,3.4067e-16,3.4069e-16,3.4071e-16,3.4073e-16,3.4075e-16,3.4076e-16,3.4078e-16,3.408e-16,3.4082e-16,3.4084e-16,3.4086e-16,3.4088e-16,3.409e-16,3.4092e-16,3.4094e-16,3.4096e-16,3.4098e-16,3.41e-16,3.4102e-16,3.4104e-16,3.4106e-16,3.4108e-16,3.411e-16,3.4112e-16,3.4114e-16,3.4116e-16,3.4119e-16,3.4121e-16,3.4123e-16,3.4125e-16,3.4127e-16,3.4129e-16,3.4131e-16,3.4133e-16,3.4135e-16,3.4137e-16,3.4139e-16,3.4141e-16,3.4143e-16,3.4145e-16,3.4147e-16,3.415e-16,3.4152e-16,3.4154e-16,3.4156e-16,3.4158e-16,3.416e-16,3.4162e-16,3.4164e-16,3.4166e-16,3.4169e-16,3.4171e-16,3.4173e-16,3.4175e-16,3.4177e-16,3.4179e-16,3.4181e-16,3.4184e-16,3.4186e-16,3.4188e-16,3.419e-16,3.4192e-16,3.4194e-16,3.4197e-16,3.4199e-16,3.4201e-16,3.4203e-16,3.4205e-16,3.4207e-16,3.421e-16,3.4212e-16,3.4214e-16,3.4216e-16,3.4218e-16,3.4221e-16,3.4223e-16,3.4225e-16,3.4227e-16,3.4229e-16,3.4231e-16,3.4234e-16,3.4236e-16,3.4238e-16,3.424e-16,3.4242e-16,3.4245e-16,3.4247e-16,3.4249e-16,3.4251e-16,3.4253e-16,3.4255e-16,3.4257e-16,3.426e-16,3.4262e-16,3.4264e-16,3.4266e-16,3.4268e-16,3.427e-16,3.4272e-16,3.4274e-16,3.4276e-16,3.4278e-16,3.428e-16,3.4282e-16,3.4284e-16,3.4286e-16,3.4288e-16,3.429e-16,3.4292e-16,3.4294e-16,3.4296e-16,3.4298e-16,3.4299e-16,3.4301e-16,3.4303e-16,3.4305e-16,3.4306e-16,3.4308e-16,3.4309e-16,3.4311e-16,3.4312e-16,3.4313e-16,3.4315e-16,3.4316e-16,3.4317e-16,3.4318e-16,3.4319e-16,3.4319e-16,3.432e-16,3.4321e-16,3.4321e-16,3.4321e-16,3.4322e-16,3.4322e-16,3.4322e-16,3.4321e-16,3.4321e-16,3.432e-16,3.4319e-16,3.4318e-16,3.4317e-16,3.4316e-16,3.4314e-16,3.4312e-16,3.431e-16,3.4307e-16,3.4304e-16,3.4301e-16,3.4298e-16,3.4294e-16,3.4289e-16,3.4285e-16,3.428e-16,3.4274e-16,3.4268e-16,3.4261e-16,3.4254e-16,3.4246e-16,3.4238e-16,3.4229e-16,3.4219e-16,3.4209e-16,3.4197e-16,3.4185e-16,3.4172e-16,3.4158e-16,3.4143e-16,3.4128e-16,3.4111e-16,3.4104e-16,3.41e-16,3.4097e-16,3.4093e-16,3.4089e-16,3.4085e-16,3.4081e-16,3.4076e-16,3.4072e-16,3.4068e-16,3.4063e-16,3.4058e-16,3.4054e-16,3.4049e-16,3.4044e-16,3.4039e-16,3.4033e-16,3.4028e-16,3.4023e-16,3.4017e-16,3.4011e-16,3.4005e-16,3.4e-16,3.3993e-16,3.3987e-16,3.3981e-16,3.3974e-16,3.3968e-16,3.3961e-16,3.3954e-16,3.3947e-16,3.394e-16,3.3933e-16,3.3925e-16,3.3918e-16,3.391e-16,3.3902e-16,3.3894e-16,3.3886e-16,3.3878e-16,3.3869e-16,3.3861e-16,3.3852e-16,3.3843e-16,3.3834e-16,3.3824e-16,3.3801e-16,3.3777e-16,3.3753e-16,3.3727e-16,3.3701e-16,3.3673e-16,3.3645e-16,3.3616e-16,3.3585e-16,3.3554e-16,3.3522e-16,3.3488e-16,3.3453e-16,3.3418e-16,3.3381e-16,3.3342e-16,3.3303e-16,3.3262e-16,3.322e-16,3.3176e-16,3.3131e-16,3.3084e-16,3.3036e-16,3.2987e-16,3.2936e-16,3.2883e-16,3.2828e-16,3.2772e-16,3.2714e-16,3.2654e-16,3.2592e-16,3.2528e-16,3.2462e-16,3.2394e-16,3.2324e-16,3.2252e-16,3.2178e-16,3.2101e-16,3.2022e-16,3.1941e-16,3.1857e-16,3.1797e-16,3.1752e-16,3.1705e-16,3.1658e-16,3.161e-16,3.1562e-16,3.1513e-16,3.1463e-16,3.1412e-16,3.136e-16,3.1308e-16,3.1255e-16,3.1201e-16,3.1146e-16,3.1091e-16,3.1035e-16,3.0977e-16,3.0919e-16,3.086e-16,3.08e-16,3.074e-16,3.0678e-16,3.0616e-16,3.0552e-16,3.0488e-16,3.0423e-16,3.0357e-16,3.0289e-16,3.0221e-16,3.0152e-16,3.0082e-16,3.0011e-16,2.9939e-16,2.9866e-16,2.9792e-16,2.9717e-16,2.964e-16,2.9563e-16,2.9536e-16,2.9516e-16,2.9497e-16,2.9478e-16,2.9459e-16,2.9439e-16,2.942e-16,2.94e-16,2.9381e-16,2.9361e-16,2.9342e-16,2.9322e-16,2.9303e-16,2.9283e-16,2.9263e-16,2.9243e-16,2.9224e-16,2.9204e-16,2.9184e-16,2.9164e-16,2.9144e-16,2.9124e-16,2.9104e-16,2.9084e-16,2.9064e-16,2.9043e-16,2.9023e-16,2.9003e-16,2.8982e-16,2.8962e-16,2.8942e-16,2.8921e-16,2.8901e-16,2.8881e-16,2.8857e-16,2.8802e-16,2.8747e-16,2.8691e-16,2.8635e-16,2.8578e-16,2.8521e-16,2.8464e-16,2.8406e-16,2.8348e-16,2.8289e-16,2.8229e-16,2.817e-16,2.811e-16,2.8049e-16,2.7988e-16,2.7927e-16,2.7865e-16,2.7802e-16,2.7739e-16,2.7676e-16,2.7612e-16,2.7548e-16,2.7483e-16,2.7418e-16,2.7352e-16,2.7286e-16,2.7219e-16,2.7152e-16,2.7084e-16,2.7016e-16,2.6947e-16,2.6878e-16,2.685e-16,2.6827e-16,2.6804e-16,2.678e-16,2.6757e-16,2.6733e-16,2.671e-16,2.6686e-16,2.6663e-16,2.6639e-16,2.6615e-16,2.6592e-16,2.6568e-16,2.6544e-16,2.652e-16,2.6497e-16,2.6474e-16,2.6449e-16,2.6425e-16,2.6402e-16,2.6378e-16,2.6354e-16,2.633e-16,2.6307e-16,2.6282e-16,2.6258e-16,2.6235e-16,2.6211e-16,2.6186e-16,2.6162e-16,2.5863e-16,2.5539e-16,2.5201e-16,2.4853e-16,2.4495e-16,2.4122e-16,2.3737e-16,2.3338e-16,2.2928e-16,2.2501e-16,2.2059e-16,2.1604e-16,2.1134e-16,2.0644e-16,2.0138e-16,1.9619e-16,1.9079e-16,1.8519e-16,1.7937e-16,1.7342e-16,1.6723e-16,1.6081e-16,1.5414e-16,1.4729e-16,1.4019e-16,1.3281e-16,1.2517e-16,1.1727e-16,1.1647e-16,1.1625e-16,1.1609e-16,1.1589e-16,1.1571e-16,1.1548e-16,1.153e-16,1.1514e-16,1.1496e-16,1.1472e-16,1.1453e-16,1.1438e-16,1.142e-16,1.1396e-16,1.1377e-16,1.1363e-16,1.1345e-16,1.1322e-16,1.13e-16,1.1286e-16,1.1271e-16,1.1249e-16,1.1225e-16,1.121e-16,1.1195e-16,1.1176e-16,1.1153e-16,1.1135e-16,1.1118e-16,1.11e-16,1.108e-16,1.1061e-16,1.1043e-16,1.1026e-16,1.1008e-16,1.0986e-16,1.0968e-16,1.0953e-16,1.0937e-16,1.0913e-16,1.0891e-16,1.0877e-16,1.0864e-16,1.0842e-16,1.0818e-16,1.08e-16,1.0788e-16,1.0769e-16,1.0748e-16,1.0727e-16,1.0713e-16,1.0674e-16,1.0623e-16,1.0573e-16,1.0528e-16,1.0484e-16,1.0433e-16,1.0382e-16,1.0333e-16,1.029e-16,1.0245e-16,1.0199e-16,1.0142e-16,1.0097e-16,1.0053e-16,1.0003e-16,9.9537e-17,9.909e-17,9.8595e-17,9.8129e-17,9.7692e-17,9.7138e-17,9.6671e-17,9.6246e-17,9.5755e-17,9.5242e-17,9.482e-17,9.4291e-17,9.3842e-17,9.3428e-17,9.2871e-17,9.2369e-17,9.1912e-17,9.1424e-17,9.0978e-17,9.0551e-17,8.9988e-17,8.9536e-17,8.9081e-17,8.856e-17,8.8132e-17,8.7599e-17,8.7107e-17,8.6686e-17,8.6196e-17,8.5672e-17,8.5214e-17} };

    template <typename Float>
	const static Float base_density = 1.8662e+18;//km^-3

    template <typename Float>
	Float Hp = 2000e-3; // Desert value

    template <typename Float>
	Float pb = 2e6; // background value in km^-3

    template <typename Float, typename Spectrum, typename Wavelength>
	class DesertDustAerosol : public GlobalAerosolModel<Float, Spectrum, Wavelength> {
	public:

        DesertDustAerosol() = default;

        Spectrum get_absorption() const override {
            return Utils::max(std::vector<Float>(tabulatedValues<Float>[0], tabulatedValues<Float>[0] + 1000), tabulatedValues<Float>[1]);
        }

        Spectrum get_absorption(const Wavelength &wl) const override {
            Spectrum s(0.);

            for (int i = 0; i < wl.Size; i++)
                s[i] = Utils::interpolate(std::vector<Float>(tabulatedValues<Float>[0], tabulatedValues<Float>[0] + 1000), tabulatedValues<Float>[1], wl[i]);

            return s;
        }

        Spectrum get_scattering() const override {
            return Utils::max(std::vector<Float>(tabulatedValues<Float>[0], tabulatedValues<Float>[0] + 1000), tabulatedValues<Float>[2]);
        }

        Spectrum get_scattering(const Wavelength &wl) const override {
            Spectrum s(0.);

            for (int i = 0; i < wl.Size; i++)
                s[i] = Utils::interpolate(std::vector<Float>(tabulatedValues<Float>[0], tabulatedValues<Float>[0] + 1000), tabulatedValues<Float>[2], wl[i]);

            return s;
        }

		Float get_density(Float z) const override {
            return base_density<Float> * (exp(-z / Hp<Float>) + pb<Float> / base_density<Float>);
        }
	};
}
#endif //_AEROSOL_DESERT_DUST