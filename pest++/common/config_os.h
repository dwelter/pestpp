#ifndef CONFIG_OS_H_
#define CONFIG_OS_H_


#define PESTPP_VERSION "3.2.0";

#if defined(_WIN32) || defined(_WIN64)
#define OS_WIN
#define DEF_DLAMCH DLAMCH
#define DEF_DLANBPRO_SPARCE DLANBPRO_SPARCE
#elif defined( __linux__)
#define OS_LINUX
#define DEF_DLAMCH dlamch_
#define DEF_DLANBPRO_SPARCE dlanbpro_sparce_
#elif defined (__APPLE__)
#define OS_LINUX
#define DEF_DLAMCH dlamch_
#define DEF_DLANBPRO_SPARCE dlanbpro_sparce_
#endif


#endif /* CONFIG_OS_H_ */
