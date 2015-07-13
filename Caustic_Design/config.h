#ifndef CONFIG_H
#define CONFIG_H
    /*Voronoi Generation*/
    #define SITE_AMOUNT 5000
    #define EPSILON 0.5
    #define LLYOD_STEPS 20

    /*Optimal Transport*/
    #define LEVEL_MAX 1

    /*Interpolation*/
    #define MESH_AMOUNT 50

    /* debug and demo */
    #define LIVE_DEMO

    /* manual optimization implementation */
    //#define DESCENT_GRADIENT

    #define DOMAIN_WIDTH 0.5
#endif // CONFIG_H
