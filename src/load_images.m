function [images, Damage_E, Damage_I, Damage_astro] = load_images(Pattern, mneuro_I)
    
    % if Pattern == 1 

    % else
        images_dir = '../images - 71';
        image_names = {
            '5_71x71.jpg', ... %1 order: 5,p,h,9,c
            'p_71x71.jpg', ... %2
            'h_71x71.jpg', ... %3
            '9_71x71.jpg', ... %4
            'c_71x71.jpg', ... %5
            'white.jpg',   ... %6
            'j_71x71.jpg', ... %7
            '8_71x71.jpg', ... %8
            'A_71x71.jpg', ... %9
            '0_71x71.jpg', ... %10
            'u_71x71.jpg', ... %11
            'e_71x71.jpg'      %12
        };
        Damage_E = imread(fullfile(images_dir, 'Damage71x71_3col_80-60.jpg'));	
        Damage_astro = imread(fullfile(images_dir, 'Damage_oval_71x71_astro.jpg'));
    % end 
    images = {};
    for name = image_names
        image = imread(fullfile(images_dir, name{1}));
        image = rgb2gray(image);
        images{end + 1} = image;
    end
       
    Damage_E = rgb2gray(Damage_E);
    Damage_I = imresize(Damage_E, [mneuro_I mneuro_I]);
    Damage_astro = rgb2gray(Damage_astro);
    
end
