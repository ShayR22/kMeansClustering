
#ifndef VECTOR_HEADER
#define VECTOR_HEADER
	
#define NUM_OF_ELEMENTS_IN_VECTOR 3

typedef struct VECTOR
{
	//axis are x,y,z
	double axis[NUM_OF_ELEMENTS_IN_VECTOR];
}vector_t;

double vector_get_distance(vector_t *v1, vector_t *v2);
void vector_add_vector(vector_t *v, vector_t *toBeAdd);
void vector_divide_vector(vector_t *v, double divide);
int vector_is_equal(vector_t *v1, vector_t *v2);

void vector_print_vector(vector_t *vector);

#endif // !VECTOR_HEADER

