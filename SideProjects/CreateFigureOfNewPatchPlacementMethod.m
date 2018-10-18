variables.patches.type = 'park';
variables.patches.size_x = 20;
variables.patches.size_y = 20;
variables.data.max_x = 200;
variables.data.max_y = 200;
variables.patches.ueyelid = 1;
variables.patches.leyelid = 1;
variables.patches.tear = 1;
variables.patches.ocorner = 1;
variables.patches.ieyebrow = 1;
variables.patches.oeyebrow = 1;
variables.patches.skin = 1;
variables.patches.side = 'left';
result = getpatches(variables);

I = zeros(200);
for i = 1:result.num
   I(result.y(i),result.x(i):result.dx(i)) = 255;
   I(result.dy(i),result.x(i):result.dx(i)) = 255;
   I(result.y(i):result.dy(i),result.x(i)) = 255;
   I(result.y(i):result.dy(i),result.dx(i)) = 255;
end

imwrite(I, 'park_patches.bmp');