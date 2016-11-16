

#ifndef _KMCODE_H_
#define _KMCODE_H_

#include <inttypes.h>

#if __cplusplus
extern "C" {
#endif

struct kmcode32_base
{
	double  origin[3],d[3];
    int32_t binary_oper[3]; //pos_in_line,line_in_layer
    int8_t in_line_digits,line_digits,layer_digits,shift;
};

struct kmcode64_base
{
	double  origin[3],d[3];
    int64_t binary_oper[3]; //pos_in_line,line_in_layer
    int8_t in_line_digits,line_digits,layer_digits,shift;
};

int kmcode32_create_base_3d(const double Min_point[3], const double Max_point[3], kmcode32_base* Base);
int kmcode32_encode_3d(const kmcode32_base* Base, const double Points[],const int N_points, int32_t * Encoded_point_ids);
int	kmcode32_decode_3d(const kmcode32_base* Base, const int32_t Encoded_point_ids[], const int N_points, double* Points);

#if __cplusplus
}
#endif

#endif