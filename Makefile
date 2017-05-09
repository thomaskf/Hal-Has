all : HAL HAS

HAL : RAL_new.h alignment.h charSet.h matrix.h tool_box.h RAL_new.cpp alignment.cpp charSet.cpp main_BU.cpp matrix.cpp optimization.h optimization.cpp core.h core.cpp rateMatrixResult.h rateMatrixResult.cpp parameters.h parameters.cpp user_options.h user_options.cpp lbfgsb_new.h lbfgsb_new.cpp definitions.h gradient.h gradient.cpp jacobi_eigenvalue.h jacobi_eigenvalue.cpp tool_box.cpp RAL_common.h RAL_common.cpp simd.h
	g++ -o HAL RAL_new.cpp alignment.cpp charSet.cpp main_BU.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp -lpthread -O3 -msse4.2 

HAS : alignment.h charSet.h matrix.h tool_box.h alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.h optimization.cpp core.h core.cpp rateMatrixResult.h rateMatrixResult.cpp parameters.h parameters.cpp RAS.h RAS.cpp user_options.h user_options.cpp lbfgsb_new.h lbfgsb_new.cpp definitions.h gradient.h gradient.cpp jacobi_eigenvalue.h jacobi_eigenvalue.cpp tool_box.cpp RAL_common.h RAL_common.cpp global.h global.cpp simd.h
	g++ -o HAS alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp RAS.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp -lpthread -O3 -msse4.2

clean :
	rm -f HAL HAS
