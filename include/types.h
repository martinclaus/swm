/*
 * Define default data types
*/

#ifndef FILE_TYPE_SEEN
#define FILE_TYPE_SEEN

/* default integer kind */
#if ! defined(INT_KIND)
#define INT_KIND 4
#endif

/* kind for mask arrays */
#if ! defined(SHORT_KIND)
#define SHORT_KIND 1
#endif

/* Kind of time index */
#if ! defined(INT_ITT_KIND)
#define INT_ITT_KIND 8
#endif

/* Kind of integer used for netcdf interface */
#if ! defined(INT_NF90_KIND)
#define INT_NF90_KIND 4
#endif


/* working precision of real variables */
#if ! defined(DOUBEL_KIND)
#define DOUBLE_KIND 8
#endif

#endif
