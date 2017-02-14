fileId = fopen('norm');
arr = fscanf(fileId, '%f', [256, 256]);
fclose(fileId);
arr = arr';

figure;
contourf(arr);
colorbar;