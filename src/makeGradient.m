function COLORS = makeGradient(numofColors, mainColors, locs, showOrNot)
	% numofColors: number of colors in output
	% mainColors : gradient colors from left to right
	% locs       : position of colors
	% showOrNot  : if 1 show the gradient, else does not show the gradient.
	
	numofGrads	= size(mainColors, 1) - 1;
	cols		= zeros(numofGrads, 1);
	areaa       = [locs(1), locs(end)];
    
	for i = 1:(numofGrads + 1)
		locs(i)	= lineMapping(locs(i), areaa, [0, 1]);
	end
   
	for i = 1:numofGrads
		cols(i) = ceil(numofColors * (locs(i + 1) - locs(i)) );
   end
    
	initcolors	= zeros(sum(cols), 3);
	shifft = 0;
	
	%gradient generator
	for i = 1:numofGrads
		color1		= mainColors(i, :);
		color2		= mainColors(i + 1, :);
		gradian		= interp1([0, 1], [color1; color2], linspace(0, 1, cols(i)));
		initcolors((1:cols(i)) + (i - 1) + shifft, :) = gradian;
		shifft		= shifft + cols(i) - 1;
	end
	COLORS = initcolors(1:numofColors, :);
	
	% Plot gradient
	if showOrNot == 1
		figure();
		img	= repmat(1:numofColors, numofColors, 1);
		img = img(:,end:-1:1)';
		imshow(img, COLORS);
		colorbar;
	end

end

function output = lineMapping(x, from, to)
	% Mapp 'x' from band [w1, v1] to band [w2, v2]
	w1 = from(1);
	v1 = from(2);
	w2 = to(1);
	v2 = to(2);
	output = 2*((x - w1)/(v1 - w1)) - 1;
	output = (output + 1)*(v2 - w2)/2 + w2;
end
