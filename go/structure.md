# adjust

## class

1. class Daoxian

    attribute:
    - input: edges, corners, and error info
    - output: error, accuracy, parameters info
    - Intermediate: init_params
    
    function:
    - check_info(): check input is structured or not. Main problems is a \n or empty line or empty but structured line.
    - init_params(): init these params. May calculate without adjustment, and get similar params?
    - cal(): main function to adjust.
    - plt_ellipsoid(): output a svg.
    - out_code(): output matlab code?
    
2. class WebApp

    attribute:
    - 
    
    function:
    - get_input():
    
    