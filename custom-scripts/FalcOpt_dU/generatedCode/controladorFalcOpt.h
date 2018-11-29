/* Copyright (c) 2017, ETH Zurich, Automatic Control Laboratory */
/* Main developers: */
/* Giampaolo Torrisi <giampaolo.torrisi@gmail.com> */
/* Damian Frick <falcopt@damianfrick.com> */
/* Tommaso Robbiani <tommasro@student.ethz.ch> */
/* Scientic Contributors: */
/* Sergio Grammatico */
/* Roy S. Smith */
/* Manfred Morari */

#ifndef _CONTROLADORFALCOPT_H_
#define _CONTROLADORFALCOPT_H_

/** FALCOPT parameters **/
#define CONTROLADOR_FALC_OPT_MAXIMUM_ITERATIONS (8000) /* Maximum number of iterations */
#define CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE (0.001) /* KKT optimality tolerance */
/* Derived parameters DO NOT TOUCH */
#define CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE_SQ (CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE*CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE)
/** Line Search Parameters **/
#define CONTROLADOR_FALC_OPT_MAXIMUM_LINE_SEARCH_ITERATIONS (50) /* Maximum number of line search iterations */
#define CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_MAX (0.4) /* Alpha max - DO NOT CHANGE, switch to variable step size to make it adjustable */
/* Derived quantities DO NOT TOUCH */
#define CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_TOL (CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_MAX*CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE) /* alphaMax*KKTTol */
#define CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_TOL_SQ (CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_TOL*CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE) /* alphaMax*KKTTol^2 */
#define CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_SQ_TOL (CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_MAX*CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_TOL) /* alphaMax^2*KKTTol */
#define CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_SQ_TOL_SQ (CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_SQ_TOL*CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE) /* alphaMax^2*KKTTol^2 */

#endif /* _CONTROLADORFALCOPT_H_ */
