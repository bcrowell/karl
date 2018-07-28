# return codes for Runge-Kutta, designed to be bitwise or-able.

#define RK_ERR 1
# ... something went really wrong, output is garbage
#define RK_INCOMPLETE 2
# ... the geodesic was incomplete
#define RK_TRIGGER 4
# ... exited due to a trigger
