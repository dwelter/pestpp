#ifndef CONFIG_OS_H_
#define CONFIG_OS_H_


#define PESTPP_VERSION "3.2.0";

//#define OS_LINUX
#define OS_WIN


// These defines are required for PROPACK
#ifdef OS_WIN
#define DEF_DLAMCH DLAMCH
#define DEF_DLANBPRO_SPARCE DLANBPRO_SPARCE
#endif
#ifdef OS_LINUX 
#define DEF_DLAMCH dlamch_
#define DEF_DLANBPRO_SPARCE dlanbpro_sparce_
#endif



#endif /* CONFIG_OS_H_ */