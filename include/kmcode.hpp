#ifndef _KMCODE_H_
#define _KMCODE_H_

#include <cinttypes>

template<typename T_PointID, int TDimension>
class KMCode;
{
	double  origin[TDimension],d[TDimension];
    T_PointID binary_oper[TDimension]; //pos_in_line,line_in_layer
    int8_t in_line_digits,line_digits,layer_digits,shift;
	
public:
	void createBase(const double Min_point[TDimension],const double Max_point[TDimension]);
	{
		const double span[3] = { max[0] - min[0], max[1] - min[1], max[2] - min[2] } ;
		const int smallest_dim = std::min_element(span,span+3) - span;
//        const int biggest_dim = std::max_element(span, span+3) - span;
		int ratio[3] = { int(span[0]/span[smallest_dim]), int(span[1]/span[smallest_dim]), int(span[2]/span[smallest_dim]) };
		int ratio_sum = ratio[0]+ratio[1]+ratio[2];

		mf_log_info("Compressed mesh ratios %d : %d : %d",ratio[0],ratio[1],ratio[2]);

		// NOTE: if failing with accuracy in smallest dim increase min_digits
		const int min_digits=6;
		int n_bits = std::numeric_limits<PTID>::digits-(3*min_digits);

//        // NOTE: we assume, that in each direction at least 5-level adaptaption is posiible (2^5)
//        const int min_n_adapts_in_dim = 5;
//        n_bits-=min_n_adapts_in_dim*3;
//        ratio[smallest_dim]=min_n_adapts_in_dim;
//        n_bits-=min_n_adapts_in_dim;

		double avg = double(n_bits)/double(ratio_sum);

		mf_log_info("Avg= %lf",avg);

		for(int i=0; i<3 ;++i) {
			ratio[i] = min_digits + double(ratio[i])*avg;
//                    + (min_n_adapts_in_dim-1);
		}

		ratio_sum = ratio[0]+ratio[1]+ratio[2];
		assert( ratio_sum <= std::numeric_limits<PTID>::digits );

		while(ratio_sum < std::numeric_limits<PTID>::digits) {
			++ratio[smallest_dim];
			++ratio_sum;
		}

		for(int i=0; i < 3; ++i) {
			origin[i] = min[i];
			d[i] = span[i] / ((1<<ratio[i])-1);
		}

		setLengths(ratio[0],ratio[1]);


		mf_log_info("Min= %lf, %lf, %lf",min[0],min[1],min[2]);
		mf_log_info("Max= %lf, %lf, %lf",max[0],max[1],max[2]);
		mf_log_info("Span= %lf, %lf, %lf",span[0],span[1],span[2]);
		mf_log_info("Smallest dim=%d",smallest_dim);
		mf_log_info("Compressed mesh ratios %d : %d : %d",ratio[0],ratio[1],ratio[2]);
		mf_log_info("pt per line=%d, lines=%d, layers=%d",in_line_digits,line_digits,layer_digits);
		mf_log_info("binary_oper[3] = %d, %d, %d",binary_oper[0],binary_oper[1],binary_oper[2]);
		mf_log_info("d[3]=%lf, %lf, %lf",d[0], d[1], d[2]);

		mf_check(min[0] == origin[0], "Error in encoding processor!");
		mf_check(min[1] == origin[1], "Error in encoding processor!");
		mf_check(min[2] == origin[2], "Error in encoding processor!");
		mf_check(max[0] == origin[0]+d[0]*((1<<ratio[0])-1), "Error in encoding processor!");
		mf_check(max[1] == origin[1]+d[1]*((1<<ratio[1])-1), "Error in encoding processor!");
		mf_check(max[2] == origin[2]+d[2]*((1<<ratio[2])-1), "Error in encoding processor!");
	}
	
	int encode(const double Points[],const int N_points, T_PointID * Encoded_point_ids);
	{
		assert(PTID((coords[0]-origin[0]) / d[0]) < (1<<in_line_digits));
        assert(PTID((coords[1]-origin[1]) / d[2])>>in_line_digits < (1<<line_digits));
        assert(PTID((coords[2]-origin[2]) / d[2])>>(in_line_digits+line_digits) < (1<<layer_digits));

        return    PTID((coords[0]-origin[0]) / d[0]) // pos in line
                | (PTID((coords[1]-origin[1]) / d[1]) << in_line_digits) // which line
                | (PTID((coords[2]-origin[2]) / d[2]) << (line_digits+in_line_digits)); // which layer

	}
	
	int decode(const T_PointID Encoded_point_ids[], const int N_points, double* Points)
	{
		coords[0]=origin[0] + d[0]*(pt & (binary_oper[0])) ;
        coords[1]=origin[1] + d[1]*((pt & (binary_oper[1]))>>in_line_digits) ;
        coords[2]=origin[2] + d[2]*(pt >> (in_line_digits+line_digits));
	}
	
private:
    void setLengths(int lineLengthAs2Pow, int layerLengthAs2Pow)
    {
        in_line_digits = lineLengthAs2Pow;
        line_digits = layerLengthAs2Pow;
        layer_digits = digits - (lineLengthAs2Pow+layerLengthAs2Pow);
        shift = in_line_digits + line_digits;

        binary_oper[0]=0;
        binary_oper[1]=0;
        binary_oper[2]=0;

        for(int i=0;i<in_line_digits; ++i) {
            binary_oper[0]|=(1<<i);
        }

        for(int i=0;i<line_digits; ++i) {
            binary_oper[1]|=(1<<i);
        }
        binary_oper[1] = binary_oper[1] << lineLengthAs2Pow;

        for(int i=0;i<layer_digits; ++i) {
            binary_oper[2]|=(1<<i);
        }
        binary_oper[2] = (binary_oper[2] << (lineLengthAs2Pow+layerLengthAs2Pow));
    }
	
};



typedef KMCode<int32_t,1>	KMCode32_1D;
typedef KMCode<int32_t,2>	KMCode32_2D;
typedef KMCode<int32_t,3>	KMCode32_3D;

typedef KMCode<int64_t,1>	KMCode64_1D;
typedef KMCode<int64_t,2>	KMCode64_2D;
typedef KMCode<int64_t,3>	KMCode64_3D;


#endif