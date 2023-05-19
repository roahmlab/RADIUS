%Class used to build pdfs given some bounding box such that the pdfgives
%the distribution of the COM of a vehicle over the bounidng box size

classdef obs_pdf
   properties 
       pdf_type = 'gaussian';
       mu = [];
       sigma = [];
       
       rot_angle = 0;
       R = @(th)[cos(th),-sin(th);sin(th) cos(th)];
       obs_box = [0;0]; % obs_box is assumed to be 2-by-n in LOCAL FRAME!
   end
   
   methods
       function obj = obs_pdf(pdf_type,varargin)
           obj.pdf_type = pdf_type;
           if nargin>1
               obj = parse_args(obj,varargin{:}) ;
           end
           switch pdf_type
               case 'gaussian'
                   if isempty(obj.mu)
                       obj.mu = obj.R(obj.rot_angle)*0.5*[max(obj.obs_box(1,:)) + min(obj.obs_box(1,:)); max(obj.obs_box(2,:)) + min(obj.obs_box(2,:))];
                   end
                   if isempty(obj.sigma)
                       L = (max(obj.obs_box(1,:)) - min(obj.obs_box(1,:)))/2;
                       W = (max(obj.obs_box(2,:)) - min(obj.obs_box(2,:)))/2;
%                        h_scale=4; %horizontal scale
%                        v_scale=70; %vertical scale
%                        obj.sigma = obj.R(obj.rot_angle)*diag([L/h_scale,W/v_scale].^2)*obj.R(obj.rot_angle)';
                       obj.sigma = obj.R(obj.rot_angle)*diag([L/3*2,W/3].^2)*obj.R(obj.rot_angle)';
                   else
                       obj.obs_box = [-3*sqrt(obj.sigma(1)), 3*sqrt(obj.sigma(1)); -3*sqrt(obj.sigma(4)), 3*sqrt(obj.sigma(4))];
                   end
               otherwise
                   error('not implemented');
           end
           
       end
       
       function plot_pdf_local(obj)
           %plots swept volume box to use to check if the pdf lies inside
           %that box
           switch obj.pdf_type
               case 'gaussian'
                   L = (max(obj.obs_box(1,:)) - min(obj.obs_box(1,:)))/2;
                   W = (max(obj.obs_box(2,:)) - min(obj.obs_box(2,:)))/2;
                   box = obj.R(obj.rot_angle)*[L,L,-L,-L,L;W,-W,-W,W,W]+obj.mu; 
                   dbox = (max(box,[],2) - min(box,[],2))/2; 
                   plot2DGaussianBlurs(obj.mu,{obj.sigma},[min(box(1,:))-dbox(1), max(box(1,:))+dbox(1); min(box(2,:))-dbox(2), max(box(2,:))+dbox(2)],0.99,200); 
                   hold on;
                   plot(box(1,:),box(2,:),'m');
               otherwise
                   error('not implemented');
           end
       end

       function h = plot_pdf_environment(obj, num_sample) % used to plot in simulator
           switch obj.pdf_type
               case 'gaussian'
                   L = (max(obj.obs_box(1,:)) - min(obj.obs_box(1,:)))/2;
                   W = (max(obj.obs_box(2,:)) - min(obj.obs_box(2,:)))/2;
                   box = obj.R(obj.rot_angle)*[L,L,-L,-L,L;W,-W,-W,W,W]+obj.mu; 
                   dbox = 1* (max(box,[],2) - min(box,[],2))/2; 
                   h = plot2DGaussianBlurs(obj.mu,{obj.sigma},[min(box(1,:))-dbox(1), max(box(1,:))+dbox(1); min(box(2,:))-dbox(2), max(box(2,:))+dbox(2)],0.99,num_sample, 1);
                   load mycolormap2.mat
                   colormap(mycolormap)
% %                    plot(obj.obs_box(2,:),obj.obs_box(1,:),'m')
               otherwise
                   error('not implemented');
           end
       end

   end
end