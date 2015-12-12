function [] = test_adjoint_1d()

% waveletfamilies('n')

extmode = 'zpd'
%extmode = 'sym'
%extmode = 'ppd'

levels = 2;

wname  = 'db1'
dwname = 'db1';
%wname  = 'db2'
%dwname = 'db2';
%wname  = 'db3'
%dwname = 'db3';
%wname  = 'db4'
%dwname = 'db4';



%wname  = 'bior4.4'
%dwname = 'rbio4.4';

%wname  = 'bior3.5'
%dwname = 'rbio3.5';

%wname  = 'rbio3.5'
%dwname = 'bior3.5';

%wname  = 'sym3'
%dwname = 'sym3';

[Lo_D, Hi_D] = wfilters(wname, 'd'); % decomp filters
lf = length(Lo_D);

%for lx=lf:25
for lx=9
   W = analysis_mat_1d(lx, wname, extmode, levels);
   %Wadj = analysis_mat_adjoint_1d(lx, wname, dwname, extmode);
   %fprintf(1, '\\|W^T-Wadj\\|_inf = %f\n', norm(W.'-Wadj,'inf'))

   W
end

end
