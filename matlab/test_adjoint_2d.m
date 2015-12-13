function [] = test_adjoint_2d()

% waveletfamilies('n')

dwtmode('zpd', 'nodisp');

extmode = 'zpd'
extmode = 'sym'
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

min_sig_len = 2^levels*(lf-1);
for lx=min_sig_len:min_sig_len+10
%for lx=3
   lX = [lx lx];
   W = analysis_mat_2d(lX, wname, extmode, levels);
   %LW = build_wavedec_levels_2d(lX, levels, wname, extmode)
   Wadj = analysis_mat_adjoint_2d(lX, wname, dwname, extmode, levels);
   
   %W.'
   %Wadj
   fprintf(1, '\\|W^T-Wadj\\|_inf = %e\n', norm(W.'-Wadj,'inf'))

   x = randn(lX);
   xe = wextend('2D', extmode, x, lf-1, 'b');
   dwtmode('zpd', 'nodisp');
   [C,L] = wavedec2(xe, levels, wname);
   fprintf(1, '\\|W*x-wavedec(xe)\\|_2 = %e\n', norm(W*x(:)-C(:)));
   
   fprintf(1, '\n');
end

end
