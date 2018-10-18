function result = getpatches(variables)
%GETPATCHES Summary of this function goes here
%   Detailed explanation goes here

	%% non overlap
	if strcmp(variables.patches.type,'nonover')

		tempx = 1:variables.patches.size_x:(variables.data.max_x-variables.patches.size_x+1);
		tempy = 1:variables.patches.size_y:(variables.data.max_y-variables.patches.size_y+1);
		offx = floor((variables.data.max_x - (tempx(length(tempx)) + variables.patches.size_x-1))/2);
		offy = floor((variables.data.max_y - (tempy(length(tempy)) + variables.patches.size_y-1))/2);
		tempx = tempx + offx; ntx = numel(tempx);
		tempy = tempy + offy; nty = numel(tempy);
		result.x = repmat(tempx,1,nty);
		result.y = reshape(repmat(tempy',1,ntx)',1,[]);
		result.dx = result.x + (variables.patches.size_x-1);
		result.dy = result.y + (variables.patches.size_y-1);
		
		result.num = size(result.x,2);
		result.size_x = variables.patches.size_x;
		result.size_y = variables.patches.size_y;
		result.xnum = length(tempx);
		result.ynum = length(tempy);
		
	%% multi-scale top down approach
	elseif strcmp(variables.patches.type,'MStop')
	
		% full image
		result.x  = [1];
		result.dx = [variables.data.max_x];
		result.y  = [1];
		result.dy = [variables.data.max_y];
		
		% 2x2 image
		tempx = [1 floor(variables.data.max_x/2)];
		tempy = [1 floor(variables.data.max_y/2)];
		offx = floor((variables.data.max_x - (floor(variables.data.max_x/2)*2))/2);
		offy = floor((variables.data.max_y - (floor(variables.data.max_y/2)*2))/2);
		tempx = tempx + offx; ntx = numel(tempx);
		tempy = tempy + offy; nty = numel(tempy);
		tempx = repmat(tempx,1,nty);
		tempy = reshape(repmat(tempy',1,ntx)',1,[]);
		result.x  = [result.x  tempx];
		result.dx = [result.dx tempx+floor(variables.data.max_x/2)-1];
		result.y  = [result.y  tempy];
		result.dy = [result.dy tempy+floor(variables.data.max_y/2)-1];
		
		% 4x4 image
		tempx = 1:floor(variables.data.max_x/4):floor(variables.data.max_x/4)*3+1;
		tempy = 1:floor(variables.data.max_y/4):floor(variables.data.max_y/4)*3+1;
		offx = floor((variables.data.max_x - (floor(variables.data.max_x/4)*4))/2);
		offy = floor((variables.data.max_y - (floor(variables.data.max_y/4)*4))/2);
		tempx = tempx + offx; ntx = numel(tempx);
		tempy = tempy + offy; nty = numel(tempy);
		tempx = repmat(tempx,1,nty);
		tempy = reshape(repmat(tempy',1,ntx)',1,[]);
		result.x  = [result.x  tempx];
		result.dx = [result.dx tempx+floor(variables.data.max_x/4)-1];
		result.y  = [result.y  tempy];
		result.dy = [result.dy tempy+floor(variables.data.max_y/4)-1];
		
		% 6x6 image
		tempx = 1:floor(variables.data.max_x/6):floor(variables.data.max_x/6)*5+1;
		tempy = 1:floor(variables.data.max_y/6):floor(variables.data.max_y/6)*5+1;
		offx = floor((variables.data.max_x - (floor(variables.data.max_x/6)*6))/2);
		offy = floor((variables.data.max_y - (floor(variables.data.max_y/6)*6))/2);
		tempx = tempx + offx; ntx = numel(tempx);
		tempy = tempy + offy; nty = numel(tempy);
		tempx = repmat(tempx,1,nty);
		tempy = reshape(repmat(tempy',1,ntx)',1,[]);
		result.x  = [result.x  tempx];
		result.dx = [result.dx tempx+floor(variables.data.max_x/6)-1];
		result.y  = [result.y  tempy];
		result.dy = [result.dy tempy+floor(variables.data.max_y/6)-1];
		
		% 8x8 image
		tempx = 1:floor(variables.data.max_x/8):floor(variables.data.max_x/8)*7+1;
		tempy = 1:floor(variables.data.max_y/8):floor(variables.data.max_y/8)*7+1;
		offx = floor((variables.data.max_x - (floor(variables.data.max_x/8)*8))/2);
		offy = floor((variables.data.max_y - (floor(variables.data.max_y/8)*8))/2);
		tempx = tempx + offx; ntx = numel(tempx);
		tempy = tempy + offy; nty = numel(tempy);
		tempx = repmat(tempx,1,nty);
		tempy = reshape(repmat(tempy',1,ntx)',1,[]);
		result.x  = [result.x  tempx];
		result.dx = [result.dx tempx+floor(variables.data.max_x/8)-1];
		result.y  = [result.y  tempy];
		result.dy = [result.dy tempy+floor(variables.data.max_y/8)-1];
		
		% 10x10 image
		tempx = 1:floor(variables.data.max_x/10):floor(variables.data.max_x/10)*9+1;
		tempy = 1:floor(variables.data.max_y/10):floor(variables.data.max_y/10)*9+1;
		offx = floor((variables.data.max_x - (floor(variables.data.max_x/10)*10))/2);
		offy = floor((variables.data.max_y - (floor(variables.data.max_y/10)*10))/2);
		tempx = tempx + offx; ntx = numel(tempx);
		tempy = tempy + offy; nty = numel(tempy);
		tempx = repmat(tempx,1,nty);
		tempy = reshape(repmat(tempy',1,ntx)',1,[]);
		result.x  = [result.x  tempx];
		result.dx = [result.dx tempx+floor(variables.data.max_x/10)-1];
		result.y  = [result.y  tempy];
		result.dy = [result.dy tempy+floor(variables.data.max_y/10)-1];
		
		% 16x16 image
		% tempx = 1:floor(variables.data.max_x/16):floor(variables.data.max_x/16)*15+1;
		% tempy = 1:floor(variables.data.max_y/16):floor(variables.data.max_y/16)*15+1;
		% offx = floor((variables.data.max_x - (floor(variables.data.max_x/16)*16))/2);
		% offy = floor((variables.data.max_y - (floor(variables.data.max_y/16)*16))/2);
		% tempx = tempx + offx; ntx = numel(tempx);
		% tempy = tempy + offy; nty = numel(tempy);
		% tempx = repmat(tempx,1,nty);
		% tempy = reshape(repmat(tempy',1,ntx)',1,[]);
		% result.x  = [result.x  tempx];
		% result.dx = [result.dx tempx+floor(variables.data.max_x/16)-1];
		% result.y  = [result.y  tempy];
		% result.dy = [result.dy tempy+floor(variables.data.max_y/16)-1];
		
		result.num = size(result.x,2);
		
	%% overlap
	elseif strcmp(variables.patches.type,'over')

		result.size_x = variables.patches.size_x;
		result.size_y = variables.patches.size_y;
		result.ynum = length(1:variables.patches.size_y:variables.data.max_y-5);
		result.xnum = length(1:variables.patches.size_x:variables.data.max_x-5);
		
		result.x = [1:variables.patches.size_x:variables.data.max_x-5 round(variables.patches.size_x/2):variables.patches.size_x:variables.data.max_x-5-round(variables.patches.size_x/2)];
		result.dx = result.x+variables.patches.size_x-1;
		result.dx(size(result.x,2)) = variables.data.max_x;
		result.x = repmat(result.x,1,result.ynum);
		result.dx = repmat(result.dx,1,result.ynum);
		
		result.y = [1:variables.patches.size_y:variables.data.max_y-5 round(variables.patches.size_y/2):variables.patches.size_y:variables.data.max_y-5-round(variables.patches.size_y/2)];
		result.dy = result.y+variables.patches.size_y-1;
		result.dy(size(result.y,2)) = variables.data.max_y;
		result.y = reshape(repmat(result.y',1,result.xnum)',1,[]);
		result.dy = reshape(repmat(result.dy',1,result.xnum)',1,[]);
		
		result.num = size(result.x,2);
	
	%% park approach
	elseif strcmp(variables.patches.type,'park')
		%0.1 patch width
		center_x = 0.5*variables.data.max_x;
		center_y = 0.5*variables.data.max_y;
		result.size_x = round(0.1*variables.data.max_x);
		result.size_y = round(0.1*variables.data.max_y);
		result.ynum = 5;
		result.xnum = 7;
		
		result.x = repmat((round(center_x-(0.05*variables.data.max_x))-(3*result.size_x)):result.size_x:round(0.75*variables.data.max_x)+round(0.05*variables.data.max_x),1,result.ynum);
		result.dx = result.x+round(0.1*variables.data.max_x)-1;
		
		result.y = (round(center_y-(0.05*variables.data.max_y))-(2*result.size_y)):result.size_x:round(0.65*variables.data.max_y)+round(0.05*variables.data.max_y);
		result.y = repmat(result.y',1,result.xnum)';
		result.y = result.y(:)';
		result.dy = result.y+round(0.1*variables.data.max_y)-1;
		
		result.num = size(result.x,2);
		
	%% miller approach
	elseif strcmp(variables.patches.type,'miller')
		
		result.size_x = -1;
		result.size_y = -1;
		% upper eyelid
        m1y=[];m1dy=[];m1x=[];m1dx=[];
		if ~isempty(variables.patches.ueyelid) && variables.patches.ueyelid == 1
			pts = pointsOnACircle([(0.25*variables.data.max_x),(0.535*variables.data.max_y)],[(0.5*variables.data.max_x),(0.425*variables.data.max_y)],[(0.75*variables.data.max_x),(0.535*variables.data.max_y)],11);
			m1y = round(pts(2,:)-round(0.025*variables.data.max_y));
            m1dy = m1y+round(0.045*variables.data.max_y);
            m1x = round(pts(1,:)-round(0.025*variables.data.max_x));
            m1dx = m1x+round(0.045*variables.data.max_x);
        end
		% lower eyelid
		m2y=[];m2dy=[];m2x=[];m2dx=[];
		if ~isempty(variables.patches.leyelid) && variables.patches.leyelid == 1
			pts = pointsOnACircle([(0.25*variables.data.max_x),(0.535*variables.data.max_y)],[(0.5*variables.data.max_x),(0.575*variables.data.max_y)],[(0.75*variables.data.max_x),(0.535*variables.data.max_y)],11);
			m2y = round(pts(2,:)-round(0.025*variables.data.max_y));
            m2dy = m2y+round(0.045*variables.data.max_y);
            m2x = round(pts(1,:)-round(0.025*variables.data.max_x));
            m2dx = m2x+round(0.045*variables.data.max_x);
        end
		% tear duct
        m3y=[];m3dy=[];m3x=[];m3dx=[];
		if ~isempty(variables.patches.tear) && variables.patches.tear == 1
			m3y = [round(0.46*variables.data.max_y) round(0.46*variables.data.max_y) round(0.535*variables.data.max_y) round(0.535*variables.data.max_y)];
            m3dy = m3y+round(0.07*variables.data.max_y);
            m3x = [round(0.175*variables.data.max_x) round(0.25*variables.data.max_x) round(0.175*variables.data.max_x) round(0.25*variables.data.max_x)];
            m3dx = m3x+round(0.07*variables.data.max_x);
			if strcmp(variables.patches.side, 'right')
				m3x = m3x + round(0.5*variables.data.max_x);
				m3dx = m3dx + round(0.5*variables.data.max_x);
			end
        end
		% outer eye corner
        m4y=[];m4dy=[];m4x=[];m4dx=[];
		if ~isempty(variables.patches.ocorner) && variables.patches.ocorner == 1
			m4y = [round(0.41*variables.data.max_y) round(0.41*variables.data.max_y) round(0.535*variables.data.max_y) round(0.535*variables.data.max_y)];
            m4dy = m4y+round(0.12*variables.data.max_y);
            m4x = [round(0.625*variables.data.max_x) round(0.75*variables.data.max_x) round(0.625*variables.data.max_x) round(0.75*variables.data.max_x)];
            m4dx = m4x+round(0.12*variables.data.max_x);
			if strcmp(variables.patches.side, 'right')
				m4x = m4x - round(0.5*variables.data.max_x);
				m4dx = m4dx - round(0.5*variables.data.max_x);
			end
        end
		% eyebrow
        m5y=[];m5dy=[];m5x=[];m5dx=[];
		m6y=[];m6dy=[];m6x=[];m6dx=[];
		if (~isempty(variables.patches.ieyebrow) && variables.patches.ieyebrow == 1) ||(~isempty(variables.patches.oeyebrow) && variables.patches.oeyebrow == 1)
			pts = pointsOnACircle([(0.125*variables.data.max_x),(0.375*variables.data.max_y)],[(0.5*variables.data.max_x),(0.25*variables.data.max_y)],[(0.875*variables.data.max_x),(0.375*variables.data.max_y)],11);
			if strcmp(variables.patches.side, 'right')
				pts(1,:) = fliplr(pts(1,:));
			end
			% inner eyebrow
			if ~isempty(variables.patches.ieyebrow) && variables.patches.ieyebrow == 1
				m5y = [round(pts(2,1:6)) round(pts(2,1:6)-round(0.075*variables.data.max_y)) round(pts(2,1:6)-round(0.15*variables.data.max_y))];
				m5dy = m5y+round(0.07*variables.data.max_y);
				m5x = [round(pts(1,1:6)-round(0.035*variables.data.max_x)) round(pts(1,1:6)-round(0.035*variables.data.max_x)) round(pts(1,1:6)-round(0.035*variables.data.max_x))];
				m5dx = m5x+round(0.07*variables.data.max_x);
			end
			% outer eyebrow
			if ~isempty(variables.patches.oeyebrow) && variables.patches.oeyebrow == 1
				m6y = [round(pts(2,7:11)) round(pts(2,7:11)-round(0.075*variables.data.max_y)) round(pts(2,7:11)-round(0.15*variables.data.max_y))];
				m6dy = m6y+round(0.07*variables.data.max_y);
				m6x = [round(pts(1,7:11)-round(0.035*variables.data.max_x)) round(pts(1,7:11)-round(0.035*variables.data.max_x)) round(pts(1,7:11)-round(0.035*variables.data.max_x))];
				m6dx = m6x+round(0.07*variables.data.max_x);
			end
        end
		% skin under the eye
        m7y=[];m7dy=[];m7x=[];m7dx=[];
		if ~isempty(variables.patches.skin) && variables.patches.skin == 1
            m7y = [repmat(round(0.65*variables.data.max_y)+1,1,7) repmat(round(0.75*variables.data.max_y)+1,1,7) repmat(round(0.85*variables.data.max_y)+1,1,7)];
            m7dy = m7y+round(0.1*variables.data.max_y)-1;
			m7x = repmat(round(0.25*variables.data.max_x)+1:round(0.1*variables.data.max_x):round(0.85*variables.data.max_x)+1,1,3);
            m7dx = m7x+round(0.1*variables.data.max_x)-1;
        end
        
        result.y = [m1y m2y m3y m4y m5y m6y m7y];
        result.dy = [m1dy m2dy m3dy m4dy m5dy m6dy m7dy];
        result.x = [m1x m2x m3x m4x m5x m6x m7x];
        result.dx = [m1dx m2dx m3dx m4dx m5dx m6dx m7dx];
        result.ynum = size(result.y, 2);
		result.xnum = size(result.x, 2);
		result.num = size(result.x,2);
		
	%% miller approach
	elseif strcmp(variables.patches.type,'miller2')
		
		result.size_x = -1;
		result.size_y = -1;
		% upper eyelid
        m1y=[];m1dy=[];m1x=[];m1dx=[];
		if ~isempty(variables.patches.ueyelid) && variables.patches.ueyelid == 1
			pts = pointsOnACircle([(0.25*variables.data.max_x),(0.535*variables.data.max_y)],[(0.5*variables.data.max_x),(0.425*variables.data.max_y)],[(0.75*variables.data.max_x),(0.535*variables.data.max_y)],7);
			m1y = round(pts(2,:)-round(0.0375*variables.data.max_y));
            m1dy = m1y+round(0.07*variables.data.max_y);
            m1x = round(pts(1,:)-round(0.0375*variables.data.max_x));
            m1dx = m1x+round(0.07*variables.data.max_x);
        end
		% lower eyelid
		m2y=[];m2dy=[];m2x=[];m2dx=[];
		if ~isempty(variables.patches.leyelid) && variables.patches.leyelid == 1
			pts = pointsOnACircle([(0.25*variables.data.max_x),(0.535*variables.data.max_y)],[(0.5*variables.data.max_x),(0.575*variables.data.max_y)],[(0.75*variables.data.max_x),(0.535*variables.data.max_y)],7);
			m2y = round(pts(2,2:6)-round(0.0375*variables.data.max_y));
            m2dy = m2y+round(0.07*variables.data.max_y);
            m2x = round(pts(1,2:6)-round(0.0375*variables.data.max_x));
            m2dx = m2x+round(0.07*variables.data.max_x);
        end
		% tear duct
        m3y=[];m3dy=[];m3x=[];m3dx=[];
		if ~isempty(variables.patches.tear) && variables.patches.tear == 1
			m3y = [round(0.46*variables.data.max_y) round(0.46*variables.data.max_y) round(0.535*variables.data.max_y) round(0.535*variables.data.max_y)];
            m3dy = m3y+round(0.07*variables.data.max_y);
            m3x = [round(0.175*variables.data.max_x) round(0.25*variables.data.max_x) round(0.175*variables.data.max_x) round(0.25*variables.data.max_x)];
            m3dx = m3x+round(0.07*variables.data.max_x);
			if strcmp(variables.patches.side, 'right')
				m3x = m3x + round(0.5*variables.data.max_x);
				m3dx = m3dx + round(0.5*variables.data.max_x);
			end
        end
		% outer eye corner
        m4y=[];m4dy=[];m4x=[];m4dx=[];
		if ~isempty(variables.patches.ocorner) && variables.patches.ocorner == 1
			m4y = [round(0.41*variables.data.max_y) round(0.41*variables.data.max_y) round(0.535*variables.data.max_y) round(0.535*variables.data.max_y)];
            m4dy = m4y+round(0.12*variables.data.max_y);
            m4x = [round(0.625*variables.data.max_x) round(0.75*variables.data.max_x) round(0.625*variables.data.max_x) round(0.75*variables.data.max_x)];
            m4dx = m4x+round(0.12*variables.data.max_x);
			if strcmp(variables.patches.side, 'right')
				m4x = m4x - round(0.5*variables.data.max_x);
				m4dx = m4dx - round(0.5*variables.data.max_x);
			end
        end
		% eyebrow
        m5y=[];m5dy=[];m5x=[];m5dx=[];
		m6y=[];m6dy=[];m6x=[];m6dx=[];
		if (~isempty(variables.patches.ieyebrow) && variables.patches.ieyebrow == 1) ||(~isempty(variables.patches.oeyebrow) && variables.patches.oeyebrow == 1)
			pts = pointsOnACircle([(0.125*variables.data.max_x),(0.325*variables.data.max_y)],[(0.5*variables.data.max_x),(0.20*variables.data.max_y)],[(0.875*variables.data.max_x),(0.325*variables.data.max_y)],7);
			if strcmp(variables.patches.side, 'right')
				pts(1,:) = fliplr(pts(1,:));
			end
			% inner eyebrow
			if ~isempty(variables.patches.ieyebrow) && variables.patches.ieyebrow == 1
				m5y = [round(pts(2,1:4))-round(0.125*variables.data.max_y) round(pts(2,1:4))];
				m5dy = m5y+round(0.12*variables.data.max_y);
				m5x = [round(pts(1,1:4)-round(0.0625*variables.data.max_x)) round(pts(1,1:4)-round(0.0625*variables.data.max_x))];
				m5dx = m5x+round(0.12*variables.data.max_x);
			end
			% outer eyebrow
			if ~isempty(variables.patches.oeyebrow) && variables.patches.oeyebrow == 1
				m6y = [round(pts(2,5:7))-round(0.125*variables.data.max_y) round(pts(2,5:7))];
				m6dy = m6y+round(0.12*variables.data.max_y);
				m6x = [round(pts(1,5:7)-round(0.0625*variables.data.max_x)) round(pts(1,5:7)-round(0.0625*variables.data.max_x))];
				m6dx = m6x+round(0.12*variables.data.max_x);
			end
        end
		% skin under the eye
        m7y=[];m7dy=[];m7x=[];m7dx=[];
		if ~isempty(variables.patches.skin) && variables.patches.skin == 1
            m7y = [repmat(round(0.60*variables.data.max_y)+1,1,10) repmat(round(0.70*variables.data.max_y)+1,1,10) repmat(round(0.80*variables.data.max_y)+1,1,10) repmat(round(0.90*variables.data.max_y)+1,1,10)];
            m7dy = m7y+round(0.1*variables.data.max_y)-1;
			m7x = repmat(round(0.00*variables.data.max_x)+1:round(0.1*variables.data.max_x):round(0.90*variables.data.max_x)+1,1,4);
            m7dx = m7x+round(0.1*variables.data.max_x)-1;
        end
        
        result.y = [m1y m2y m3y m4y m5y m6y m7y];
        result.dy = [m1dy m2dy m3dy m4dy m5dy m6dy m7dy];
        result.x = [m1x m2x m3x m4x m5x m6x m7x];
        result.dx = [m1dx m2dx m3dx m4dx m5dx m6dx m7dx];
        result.ynum = size(result.y, 2);
		result.xnum = size(result.x, 2);
		result.num = size(result.x,2);
		
		
		
		%% miller approach
	elseif strcmp(variables.patches.type,'miller3')
		
		result.size_x = -1;
		result.size_y = -1;
		% upper eyelid
        m1y=[];m1dy=[];m1x=[];m1dx=[];
		if ~isempty(variables.patches.ueyelid) && variables.patches.ueyelid == 1
			pts = pointsOnACircle([(0.25*variables.data.max_x),(0.535*variables.data.max_y)],[(0.5*variables.data.max_x),(0.425*variables.data.max_y)],[(0.75*variables.data.max_x),(0.535*variables.data.max_y)],7);
			m1y = round(pts(2,:)-round(0.0375*variables.data.max_y));
            m1dy = m1y+round(0.07*variables.data.max_y);
            m1x = round(pts(1,:)-round(0.0375*variables.data.max_x));
            m1dx = m1x+round(0.07*variables.data.max_x);
        end
		% lower eyelid
		m2y=[];m2dy=[];m2x=[];m2dx=[];
		if ~isempty(variables.patches.leyelid) && variables.patches.leyelid == 1
			pts = pointsOnACircle([(0.25*variables.data.max_x),(0.535*variables.data.max_y)],[(0.5*variables.data.max_x),(0.575*variables.data.max_y)],[(0.75*variables.data.max_x),(0.535*variables.data.max_y)],7);
			m2y = round(pts(2,2:6)-round(0.0375*variables.data.max_y));
            m2dy = m2y+round(0.07*variables.data.max_y);
            m2x = round(pts(1,2:6)-round(0.0375*variables.data.max_x));
            m2dx = m2x+round(0.07*variables.data.max_x);
        end
		% tear duct
        m3y=[];m3dy=[];m3x=[];m3dx=[];
		if ~isempty(variables.patches.tear) && variables.patches.tear == 1
			m3y = [round(0.46*variables.data.max_y) round(0.46*variables.data.max_y) round(0.535*variables.data.max_y) round(0.535*variables.data.max_y)];
            m3dy = m3y+round(0.07*variables.data.max_y);
            m3x = [round(0.175*variables.data.max_x) round(0.25*variables.data.max_x) round(0.175*variables.data.max_x) round(0.25*variables.data.max_x)];
            m3dx = m3x+round(0.07*variables.data.max_x);
			if strcmp(variables.patches.side, 'right')
				m3x = m3x + round(0.5*variables.data.max_x);
				m3dx = m3dx + round(0.5*variables.data.max_x);
			end
        end
		% outer eye corner
        m4y=[];m4dy=[];m4x=[];m4dx=[];
		if ~isempty(variables.patches.ocorner) && variables.patches.ocorner == 1
			m4y = [round(0.41*variables.data.max_y) round(0.41*variables.data.max_y) round(0.535*variables.data.max_y) round(0.535*variables.data.max_y)];
            m4dy = m4y+round(0.12*variables.data.max_y);
            m4x = [round(0.625*variables.data.max_x) round(0.75*variables.data.max_x) round(0.625*variables.data.max_x) round(0.75*variables.data.max_x)];
            m4dx = m4x+round(0.12*variables.data.max_x);
			if strcmp(variables.patches.side, 'right')
				m4x = m4x - round(0.5*variables.data.max_x);
				m4dx = m4dx - round(0.5*variables.data.max_x);
			end
        end
		% eyebrow
        m5y=[];m5dy=[];m5x=[];m5dx=[];
		m6y=[];m6dy=[];m6x=[];m6dx=[];
		if (~isempty(variables.patches.ieyebrow) && variables.patches.ieyebrow == 1) ||(~isempty(variables.patches.oeyebrow) && variables.patches.oeyebrow == 1)
			pts = pointsOnACircle([(0.125*variables.data.max_x),(0.325*variables.data.max_y)],[(0.5*variables.data.max_x),(0.20*variables.data.max_y)],[(0.875*variables.data.max_x),(0.325*variables.data.max_y)],7);
			if strcmp(variables.patches.side, 'right')
				pts(1,:) = fliplr(pts(1,:));
			end
			% inner eyebrow
			if ~isempty(variables.patches.ieyebrow) && variables.patches.ieyebrow == 1
				m5y = [round(pts(2,1:4))-round(0.125*variables.data.max_y) round(pts(2,1:4))];
				m5dy = m5y+round(0.12*variables.data.max_y);
				m5x = [round(pts(1,1:4)-round(0.0625*variables.data.max_x)) round(pts(1,1:4)-round(0.0625*variables.data.max_x))];
				m5dx = m5x+round(0.12*variables.data.max_x);
			end
			% outer eyebrow
			if ~isempty(variables.patches.oeyebrow) && variables.patches.oeyebrow == 1
				m6y = [round(pts(2,5:7))-round(0.125*variables.data.max_y) round(pts(2,5:7))];
				m6dy = m6y+round(0.12*variables.data.max_y);
				m6x = [round(pts(1,5:7)-round(0.0625*variables.data.max_x)) round(pts(1,5:7)-round(0.0625*variables.data.max_x))];
				m6dx = m6x+round(0.12*variables.data.max_x);
			end
        end
		% skin under the eye
        m7y=[];m7dy=[];m7x=[];m7dx=[];
		if ~isempty(variables.patches.skin) && variables.patches.skin == 1
            m7y = [repmat(round(0.65*variables.data.max_y)+1,1,7) repmat(round(0.75*variables.data.max_y)+1,1,7) repmat(round(0.85*variables.data.max_y)+1,1,7)];
            m7dy = m7y+round(0.1*variables.data.max_y)-1;
			m7x = repmat(round(0.25*variables.data.max_x)+1:round(0.1*variables.data.max_x):round(0.85*variables.data.max_x)+1,1,3);
            m7dx = m7x+round(0.1*variables.data.max_x)-1;
        end
		% typical nonover pattern
		m8y=[];m8dy=[];m8x=[];m8dx=[];
		if ~isempty(variables.patches.full) && variables.patches.full == 1
			yt = length(1:variables.patches.size_y:variables.data.max_y-5);
			xt = length(1:variables.patches.size_x:variables.data.max_x-5);
		
			m8x = 1:variables.patches.size_x:variables.data.max_x-5;
			m8dx = m8x+variables.patches.size_x-1;
			m8dx(size(m8x,2)) = variables.data.max_x;
			m8x = repmat(m8x,1,yt);
			m8dx = repmat(m8dx,1,yt);
			
			m8y = 1:variables.patches.size_y:variables.data.max_y-5;
			m8dy = m8y+variables.patches.size_y-1;
			m8dy(size(m8y,2)) = variables.data.max_y;
			m8y = repmat(m8y,1,xt);
			m8dy = repmat(m8dy,1,xt);
		end
        
        result.y = [m1y m2y m3y m4y m5y m6y m7y m8y];
        result.dy = [m1dy m2dy m3dy m4dy m5dy m6dy m7dy m8dy];
        result.x = [m1x m2x m3x m4x m5x m6x m7x m8x];
        result.dx = [m1dx m2dx m3dx m4dx m5dx m6dx m7dx m8dx];
        result.ynum = size(result.y, 2);
		result.xnum = size(result.x, 2);
		result.num = size(result.x,2);
		
	end

end
