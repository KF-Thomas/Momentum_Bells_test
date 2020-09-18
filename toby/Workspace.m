%A = [phi; slist]
data = fopen('R8DE-6B01C.txt','w');
fprintf(data,'%6s %12s\n','\phi','E(\phi)');
%fprintf(data,'%6.2f %12.8f\n',A);
fclose(data);