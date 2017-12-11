%load edgeVol generated from machine segmentation and write edge point info
%into xml files, so that we can have access to control points in iphotodraw
%for fine tuning ground truth


clear all; clc; close all;

% pathname = '/Users/kongj/Glioma/cellMigration3DScaffold/data/LN229LIP z slices'; outputpath = [pathname '_Interp_Smooth/'];
%pathname = '/Users/kongj/Glioma/cellMigration3DScaffold/data/LN229VC z slices'; outputpath = [pathname '_Interp_Smooth/'];
pathname = '.\NucleiSegLS\nu\3000\'; outputpath = [pathname '\'];

%s = 'minVol_50_globalCoef_0.6_localCoef_1';
%s = 'minVol_10_globalCoef_0.6_localCoef_1';
% s = 'minVol_10_globalCoef_0.6_localCoef_0.95';
% s = 'minVol_10_globalCoef_0.6_localCoef_0.9';
% s = 'minVol_10_globalCoef_0.6_localCoef_0.8';
% s = 'minVol_10_globalCoef_0.6_localCoef_0.7';
% s = 'minVol_10_globalCoef_0.6_localCoef_0.6';
% s = 'minVol_10_globalCoef_0.55_localCoef_1';
% s = 'minVol_10_globalCoef_0.55_localCoef_0.95';
% s = 'minVol_10_globalCoef_0.55_localCoef_0.9';
% s = 'minVol_10_globalCoef_0.55_localCoef_0.8';
% s = 'minVol_10_globalCoef_0.55_localCoef_0.7';
% s = 'minVol_10_globalCoef_0.55_localCoef_0.6';
% s = 'minVol_10_globalCoef_0.5_localCoef_1';
% s = 'minVol_10_globalCoef_0.5_localCoef_0.95';
% s = 'minVol_10_globalCoef_0.5_localCoef_0.9';
% s = 'minVol_10_globalCoef_0.5_localCoef_0.8';
% s = 'minVol_10_globalCoef_0.5_localCoef_0.7';
% s = 'minVol_10_globalCoef_0.5_localCoef_0.6';
% s = 'minVol_10_globalCoef_0.4_localCoef_1';
% s = 'minVol_10_globalCoef_0.3_localCoef_1';
% s = 'minVol_10_globalCoef_0.2_localCoef_1';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.95';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.9';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.8';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.7';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.6';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.6_ITER_1';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.7_ITER_1';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.8_ITER_1';
% 
% s = 'minVol_10_globalCoef_0.2_localCoef_0.6_ITER_0';
% s = 'minVol_10_globalCoef_0.2_localCoef_0.8_ITER_0';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_0';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.8_ITER_0';
% s = 'minVol_80_globalCoef_0.2_localCoef_0.6_ITER_0';
% s = 'minVol_80_globalCoef_0.2_localCoef_0.8_ITER_0';
% s = 'minVol_100_globalCoef_0.2_localCoef_0.6_ITER_0';
% s = 'minVol_100_globalCoef_0.2_localCoef_0.8_ITER_0';
% 
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_2';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_3';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_4';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_5';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_6';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_7';
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_8';
% 
% s = 'minVol_50_globalCoef_0.2_localCoef_0.6_ITER_2_watershed';

files = dir([pathname '\*.mat']);
imgnames = {files.name}.';
n = length(imgnames);
% load([pathname '/finalLabel_' s '.mat'],  'finalLabel', 'waterRidge');
imgIdx = [1,3,4,7,8,11,13,14,16,19,21,23,26];
% imgIdx = 1;
for i = 11:13
    s = sprintf('%02d_boundary_coordinates', imgIdx(i));
    load([pathname '\' s '.mat']);
    % NumObj = max(finalLabel(:));
%     NumObj = length(boundaryCoordinate(:));
%     colors = jet(NumObj);
%     colors = max(min(round(colors*255), 255), 0);
%     y = randsample(NumObj, NumObj);
%     colors = colors(y, :);

%     if exist([outputpath s], 'dir')== 0
%         mkdir([outputpath s]);
%     end

%     se = strel('disk',2);
    
%     for i = 1:size(boundaryCoordinate,3)
%         boundaryFile = sprintf('%02d.tif', i);
%         I = imread([outputpath boundaryFile]);
%         red = I(:,:,1);
%         green = I(:,:,2);
%         blue = I(:,:,3);

    

        xmlFile = sprintf('%s%02d_data.xml', outputpath,imgIdx(i));
        tifFile = sprintf('%02d_New.tif', imgIdx(i));
        try
            os = java.io.FileOutputStream(xmlFile);
            bstream = java.io.BufferedOutputStream(os);
        catch exception
        end
    
        header = ['<!--Document-->', char(10),...
            '      <Document FileVersion="1.0">', char(10),...
            '          <ExportImageSettings FileName="' tifFile '" JpegQuality="95" ResizeMode="KeepOriginalSize" ImageZoomFactor="100" ImageResizeWidth="0" ImageResizeHeight="0" ImageType="IMG_TIF" />',char(10),...
            '          <ImageOptions IsNegative="False" IsGrayscale="False" IsSepia="False" Rotation="0">',char(10),...
            '              <BrightnessContrast Brightness="0" Contrast="0" />',char(10),...
            '              <HueSaturation Hue="0" Saturation="0" Lightness="0" />',char(10),...
            '              <Canvas>',char(10),...
            '                  <Box Left="0" Top="0" Width="0" Height="0" />',char(10),...
            '                  <BackColor Alpha="255" R="255" G="255" B="255" />',char(10),...
            '              </Canvas>',char(10),...
            '              <Flip HorizontalFlip="False" VerticalFlip="False" />',char(10),...
            '          </ImageOptions>',char(10),...
            '          <Layers>',char(10),...
            '              <Layer Name="Layer1" Visible="True" LockedShapesIndex="">',char(10),...
            '                  <Shapes>',char(10),...
            ];
        bstream.write(uint8(header), 0, length(header));
    
    
    
    
    
%         slide = finalLabel(:,:,i);
%         N = max(slide(:));
%         for j = 1:N
% %             ridge = imdilate((slide==j), se) & waterRidge(:,:,i);
% %             B = bwboundaries(slide==j | ridge, 'noholes');
        B = boundaryCoordinate;
%             if isempty(B)
%                 fprintf('%s, label=%d has 0 boundary\n', boundaryFile, j);
%                 continue;
%             end
%             if length(B) > 1
%                 fprintf('%s, label=%d has %d boundaries\n', boundaryFile, j, length(B));
%             end

            for k = 1:length(B)
                b = B{k}{1,1};
%                 ind = sub2ind(size(green), b(:,1), b(:,2));
%                 red(ind) = colors(j,1);
%                 green(ind) = colors(j,2);
%                 blue(ind) = colors(j,3);


                pts = [];
                b = b - 1; %iPhotoDraw starts with (0, 0)


                if size(b,1) > 10 && size(b,1) <= 20
                    b = b(1:2:end, :);
                end
                if size(b,1) > 20 && size(b,1) <= 40
                    b = b(1:3:end, :);
                end
                if size(b,1) > 40
                    b = b(1:4:end, :);
                end



                for p = 1:size(b,1)
                    pts = [pts '                                      <Point X="' num2str(b(p,2)) '" Y="' num2str(b(p,1)) '" />',char(10)];
                end
                shape =  ['                        <Shape Type="SmoothPolygon">',char(10),...
                    '                              <Settings>',char(10),...
                    '                                  <MiscSettings GroupRendering="Unknown" />',char(10),...
                    '                                  <Font Name="Arial" Size="12" Style="Regular">',char(10),...
                    '                                      <Color Alpha="255" R="0" G="0" B="0" />',char(10),...
                    '                                  </Font>',char(10),...
                    '                                  <Line Width="1" Dash="Solid" DashOffset="0" Join="Round" OutlineType="Color">',char(10),...
                    '                                      <Color Alpha="255" R="0" G="255" B="0" />',char(10),...
                    '                                      <StartArrowHead Type="None" HeightFactor="1" WidthFactor="2" />',char(10),...
                    '                                      <EndArrowHead Type="None" HeightFactor="1" WidthFactor="2" />',char(10),...
                    '                                  </Line>',char(10),...
                    '                                  <Fill FillType="Color">',char(10),...
                    '                                      <Color Alpha="0" R="255" G="255" B="255" />',char(10),...
                    '                                      <GradientSettings Type="Linear" Angle="0" HorizontalOffset="0" VerticalOffset="0">',char(10),...
                    '                                          <StartingColor Alpha="255" R="0" G="0" B="0" />',char(10),...
                    '                                          <EndingColor Alpha="255" R="255" G="255" B="255" />',char(10),...
                    '                                          <Blend />',char(10),...
                    '                                      </GradientSettings>',char(10),...
                    '                                      <EmbeddedImage Align="Center" ImageFillType="Stretch" Alpha="255" FileName="">',char(10),...
                    '                                          <StretchSettings Type="KeepOriginalSize" Align="Center" ZoomFactor="100">',char(10),...
                    '                                              <Offset X="0" Y="0" />',char(10),...
                    '                                          </StretchSettings>',char(10),...
                    '                                          <TileSettings WrapMode="Tile">',char(10),...
                    '                                              <Offset X="0" Y="0" />',char(10),...
                    '                                          </TileSettings>',char(10),...
                    '                                          <ImageData><![CDATA[]]></ImageData>',char(10),...
                    '                                      </EmbeddedImage>',char(10),...
                    '                                  </Fill>',char(10),...
                    '                                  <TextEffect UseTextEffect="False" />',char(10),...
                    '                                  <EffectSettings>',char(10),...
                    '                                      <Shadow UseShadow="False" Angle="45" Offset="5" Size="100" BlurLevel="0">',char(10),...
                    '                                          <Color Alpha="255" R="0" G="0" B="0" />',char(10),...
                    '                                      </Shadow>',char(10),...
                    '                                      <Glow UseGlow="False" BlurLevel="20" Thickness="8">',char(10),...
                    '                                          <Color Alpha="255" R="29" G="199" B="244" />',char(10),...
                    '                                      </Glow>',char(10),...
                    '                                  </EffectSettings>',char(10),...
                    '                              </Settings>',char(10),...
                    '                              <BlockText Align="Center" VerticalAlign="Middle" RightToLeft="Unknown">',char(10),...
                    '                                 <Text></Text>',char(10),...
                    '                                  <Margin Left="0" Top="0" Right="0" Bottom="0" />',char(10),...
                    '                              </BlockText>',char(10),...
                    '                              <Data Rotation="0">',char(10),...
                    '                                  <Points>',char(10),...
                    pts,...
                    '                                  </Points>',char(10),...
                    '                              </Data>',char(10),...
                    '                          </Shape>',char(10),...
                    ];

                bstream.write(uint8(shape), 0, length(shape));
            end

%         end

        tail = ['                  </Shapes>',char(10),...
            '              </Layer>',char(10),...
            '           </Layers>',char(10),...
            '           <Snapshots />',char(10),...
            '</Document>',...
            ];
        bstream.write(uint8(tail), 0, length(tail));

        try
            bstream.flush();
            os.close();
        catch exception
        end
    
    
%         imwrite(cat(3, red, green, blue), [outputpath s '/' boundaryFile]);
end
    


% save([outputpath s '/colors.mat'], 'colors');