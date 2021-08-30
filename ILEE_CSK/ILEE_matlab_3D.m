function output_file = ILEE_matlab_3D (io_directory, stack_num, height, width, xy_unit, z_unit, k, grad_thres, tol, use_GPU)
%%data input
input_file = strcat(io_directory, 'temp_ILEE_input.mat');
img_pp_flt = importdata(input_file).img_pp_flt(:);
img_grad_flt = importdata(input_file).img_grad_flt(:);

height = double(height);
width = double(width);
stack_num = double(stack_num);
xy_unit = double(xy_unit);
z_unit = double(z_unit);
k = double(k);
grad_thres = double(grad_thres);
tol = double(tol);
use_GPU = logical(use_GPU);

global pL saved_img_para;
current_img_para = [stack_num, height, width, xy_unit, z_unit];
flt_size = size(img_pp_flt,1);
stack_size = height*width;
rescaling = xy_unit/z_unit;

if not(isequal(saved_img_para, current_img_para))
    %%making pL matrix
    row = zeros([flt_size*7,1]);
    col = zeros([flt_size*7,1]);
    value = zeros([flt_size*7,1]);
    count=0;
    
    for n = 1:flt_size
        sum_weight = 0;
        i = fix((n-1)/stack_size) + 1;
        in_plane_index = mod((n-1),stack_size) + 1;
        ii = fix((in_plane_index-1)/width) + 1;
        iii = mod(n-1, width) + 1;
        
        if(iii>1)
            count = count+1;
            row(count) = n;
            col(count) = stack_size*(i-1) + width*(ii-1) + (iii-1);
            value(count) = -1;
            sum_weight = sum_weight+1;
        end
        if(iii<width)
            count = count+1;
            row(count) = n;
            col(count) = stack_size*(i-1) + width*(ii-1) + (iii+1);
            value(count) = -1;
            sum_weight = sum_weight+1;
        end
        if(ii>1)
            count = count+1;
            row(count) = n;
            col(count) = stack_size*(i-1) + width*(ii-2) + iii;
            value(count) = -1;
            sum_weight = sum_weight+1;
        end
        if(ii<height)
            count = count+1;
            row(count) = n;
            col(count) = stack_size*(i-1) + width*ii + iii;
            value(count) = -1;
            sum_weight = sum_weight+1;
        end
        if(i>1)
            count = count+1;
            row(count) = n;
            col(count) = stack_size*(i-2) + width*(ii-1) + iii;
            value(count) = -rescaling;
            sum_weight = sum_weight + rescaling;
        end
        if(i<stack_num)
            count = count+1;
            row(count) = n;
            col(count) = stack_size*i + width*(ii-1) + iii;
            value(count) = -rescaling;
            sum_weight = sum_weight + rescaling;
        end
        
        count = count+1;
        row(count) = n;
        col(count)= stack_size*(i-1) + width*(ii-1) + iii;
        value(count)= sum_weight;
    end
saved_img_para = current_img_para;
pL = sparse(row(1:count), col(1:count), value(1:count), flt_size, flt_size);    
end

row = zeros([flt_size*7,1]);
col = zeros([flt_size*7,1]);
value = zeros([flt_size*7,1]);
count=0;
for n = 1:flt_size
    if(img_grad_flt(n)>=grad_thres)
        count = count+1;
        row(count) = n;
        col(count)= n;
        value(count)= 1;
    end
end

S = sparse(row(1:count), col(1:count), value(1:count), flt_size, flt_size);        
A=k*pL+S;
b=S*img_pp_flt;

if (use_GPU)
    reset(gpuDevice);
    A_gpu=gpuArray(A);
    b_gpu=gpuArray(b);
    result_flt_gpu = pcg(A_gpu, b_gpu, tol, 10000);
    result_flt = gather(result_flt_gpu);
else
    result_flt = pcg(A, b, tol, 10000);
end

output_file = strcat(io_directory, 'temp_ILEE_output.mat');
save(output_file, 'result_flt');
end
