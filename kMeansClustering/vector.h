
#ifndef VECTOR_HEADER
#define VECTOR_HEADER
	
#define NUM_OF_ELEMENTS_IN_VECTOR 3
#define DOUBLE_EPSILON_PRECISION 0.0001

typedef struct VECTOR
{
	//axis are x,y,z
	double axis[NUM_OF_ELEMENTS_IN_VECTOR];
}vector_t;

/*
find ecludian distance between two vectors

Parameters:
- v1: pointer to vector
- v2: pointer to vector

Returns:
- double: distance ecludian distance between the two vectors
*/
double vector_get_distance(vector_t *v1, vector_t *v2);

/*
add vector v2 to vector v1

Parameters:
- v1: pointer to vector
- v2: pointer to vector

Returns: void
*/
void vector_add_vector(vector_t *v, vector_t *toBeAdd);

/*
devide vector v1 values with corresponding vector v2 values

Parameters:
- v1: pointer to vector
- v2: pointer to vector

Returns: void
*/
void vector_divide_vector(vector_t *v, double divide);

/*
compare corresponding values between v1 and v2 using DOUBLE_EPSILON_PRECISION header const.
Because doubles are used, a precision degree wanted can be changed meaning:
what is considerd equal 6 numbers after the points for example.
the check is used by if |value1 - value2| < DOUBLE_EPSILON_PRECISION then the values are
considered to be equals

Parameters: 
- v1: pointer to vector
- v2: pointer to vector

Returns:
- 1: equal
- 0: not equal
*/
int vector_is_equal(vector_t *v1, vector_t *v2);

void vector_print_vector(vector_t *vector);

#endif // !VECTOR_HEADER

