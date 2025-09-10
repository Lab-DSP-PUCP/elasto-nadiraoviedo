function [grad_phi_sin, grad_phi_cos] = compute_phase_gradient(y, e_i)
%INPUTS
%y = senial analitica
% e_i = step para demodulacion
    % Asegúrate de que el vector y tenga suficientes elementos
    if length(y) < 2 * e_i + 1
        error('El vector de señal es demasiado corto para calcular el gradiente');
    end
    
    
    % Inicializar los vectores para el rango válido
    grad_phi_sin = zeros(1, length(y));
    grad_phi_cos = zeros(1, length(y));
    
    for x = e_i + 1 : length(y) - e_i
        % Gradiente de fase usando sin^{-1}
        grad_phi_sin(x - e_i) = asin((y(x + e_i) - y(x - e_i)) / (2 *1i * y(x)));
        
        % Gradiente de fase usando cos^{-1}
        grad_phi_cos(x - e_i) = acos((y(x + e_i) + y(x - e_i)) / (2 * y(x)));
    end

    
    grad_phi_sin(1)=grad_phi_sin(2);
    grad_phi_sin(end)=grad_phi_sin(end-2);
    grad_phi_sin(end-1)=grad_phi_sin(end-2);

    grad_phi_cos(1)=grad_phi_cos(2);
    grad_phi_cos(end)=grad_phi_cos(end-2);
    grad_phi_cos(end-1)=grad_phi_cos(end-2);

    windowSize = 3;  % longitud de la ventana (tiene que ser impar)
    grad_phi_sin = medfilt1(abs(grad_phi_sin), windowSize);
    grad_phi_cos = medfilt1(real(grad_phi_cos), windowSize);

end
