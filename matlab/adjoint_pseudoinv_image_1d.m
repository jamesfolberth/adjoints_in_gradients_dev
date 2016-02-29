function [] = adjoint_pseudoinv_1d()

% waveletfamilies('n')

dwtmode('zpd', 'nodisp');

extmode = 'zpd'
extmode = 'sym'
%extmode = 'ppd'

levels = 3;

%wname  = 'db1'
%dwname = 'db1';
%wname  = 'db2'
%dwname = 'db2';
%wname  = 'db3'
%dwname = 'db3';
%wname  = 'db4'
%dwname = 'db4';

%wname  = 'bior2.2'
%dwname = 'rbio2.2'; 

wname  = 'bior4.4'
dwname = 'rbio4.4';

%wname  = 'bior3.5'
%dwname = 'rbio3.5';

%wname  = 'rbio3.5'
%dwname = 'bior3.5';

%wname  = 'sym3'
%dwname = 'sym3';

[Lo_D, Hi_D] = wfilters(wname, 'd'); % decomp filters
lf = length(Lo_D);


min_sig_len = 2^levels*(lf-1);
for lx=min_sig_len
   W = analysis_mat_1d(lx, wname, extmode, levels);
   Wadj = analysis_mat_adjoint_1d(lx, wname, dwname, extmode, levels);
   
   %W.'
   %Wadj
   fprintf(1, '\\|W^T-Wadj\\|_inf = %e\n', norm(W.'-Wadj,'inf'))
   
   Id = eye(lx);
   Ir = eye(size(W,1));
   
   %diag(Wadj*W)

   imagesc(Wadj*W); colorbar() % analysis -> synthesis
   title('W^*W')
   %imagesc(W*Wt); colorbar() % synthesis -> analysis
   %title('WW^T')
   
   fprintf(1, '\\|I-W^*W\\|_2 = %f\n', norm(Id-Wadj*W))
   fprintf(1, '\\|I-WW^*\\|_2 = %f\n', norm(Ir-W*Wadj))

   %x = randn(lx,1);
   %xe = wextend('1D', extmode, x, lf-1, 'b');
   %dwtmode('zpd', 'nodisp');
   %[C,L] = wavedec(xe, levels, wname);
   %fprintf(1, '\\|W*x-wavedec(xe)\\|_2 = %e\n', norm(W*x-C));
   %
   %fprintf(1, '\n');
end

end
