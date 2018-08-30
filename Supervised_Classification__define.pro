;+
; 监督分类
; 最短距离,最大似然,马氏距离,SAM
;-
function Supervised_Classification::ReadParamfile, $
                           paramfile = paramfile, $
                           xmlfile = xmlfile, $
                           method = method, $
                           unclassifiedfile = unclassifiedfile, $
                           classifiedfile = classifiedfile, $
                           error = error
   catch, error_status
   if error_status ne 0 then begin
     error = !ERROR_STATE.MSG
     WriteLog, !ERROR_STATE.MSG, 0
     catch, /cancel
     return, 0
   endif                           
                           
  if ~file_test(paramfile) then begin
    error = paramfile + '参数文件不存在!'
    WriteLog, paramfile + '参数文件不存在!', 0
    print, paramfile + '参数文件不存在!'
    return, 0
  endif
  
  if file_lines(paramfile) lt 4 then begin
    error = paramfile + '参数文件格式异常!'
    WriteLog, paramfile + '参数文件格式异常!', 0
    print, paramfile + '参数文件格式异常!'
    return, 0
  endif
  
  strs = strarr(4)
  openr, lun, paramfile, /GET_LUN
  readf, lun, strs
  free_lun, lun
  
  if strmatch(strs[0], '*=*') then begin
    xmlfile = (strsplit(strs[0], '=', /extract))[1]
  endif else begin
    error = '参数文件格式异常!'
    Writelog, '参数文件格式异常!', 0
    print, '参数文件格式异常!'
    return, 0
  endelse  
  if strmatch(strs[1], '*=*') then begin
    method = (strsplit(strs[1], '=', /extract))[1]
  endif else begin
    error = '参数文件格式异常!'
    WriteLog, '参数文件格式异常!', 0
    print, '参数文件格式异常!'    
    return, 0
  endelse
  if strmatch(strs[2], '*=*') then begin
    unclassifiedfile = (strsplit(strs[2], '=', /extract))[1]
  endif else begin
    error = '参数文件格式异常!'
    WriteLog, '参数文件格式异常!', 0
    print, '参数文件格式异常!'   
    return, 0
  endelse
  if strmatch(strs[3], '*=*') then begin
    classifiedfile = (strsplit(strs[3], '=', /extract))[1]
  endif else begin
    error = '参数文件格式异常!'
    WriteLog, '参数文件格式异常!', 0
    print, '参数文件格式异常!'   
    return, 0
  endelse
  
  return, 1
end


function Supervised_Classification::ReadROIXML, $
                         xmlfile = xmlfile, $
                         num_class = num_class, $
                         class_names = class_names, $
                         lookup = lookup, $
                         points = points, $
                         error = error
  compile_opt idl2
  
  catch, error_status
  if error_status ne 0 then begin
    error = !ERROR_STATE.MSG
    WriteLog, !ERROR_STATE.MSG, 0
    catch, /cancel
    return, 0
  endif
  
  if ~file_test(xmlfile) then begin
    error = xmlfile+'ROI样本文件不存在!'
    writelog, xmlfile+'ROI样本文件不存在!', 0
    print, xmlfile+'ROI样本文件不存在!'
    return, 0
  endif
  
  odoc = obj_new('IDLffXMLDOMDocument', filename = xmlfile)
  oxmlele = odoc->GetDocumentElement()
  oregions = oxmlele->GetElementsByTagName('Region')
  num_regoins = oregions->GetLength()
  if num_regoins le 0 then begin
    error = '样本文件中不包含样本记录!'
    obj_destroy, [odoc, oxmlele, oregions]
    writelog, '样本文件中不包含样本记录!', 0
    print, '样本文件中不包含样本记录!'
    return, 0
  endif
  
  num_class = num_regoins
  points = ptrarr(num_class)
  class_names = strarr(num_regoins)
  lookup = lonarr(3, num_class)
  
  for i=0,num_regoins-1 do begin
    oregion = oregions->item(i)
    class_names[i] = oregion->GetAttribute('name')
    lookup[*,i] = long(strsplit(oregion->GetAttribute('color'), ',', /extract))
    
    ogeometrydefs = oregion->GetElementsByTagName('GeometryDef')
    if ogeometrydefs->GetLength() ne 1 then begin
      obj_destroy, [odoc, oxmlele, oregions, oregion, ogeometrydefs]
      error = '样本文件GeometryDef节点格式异常!'
      writelog, '样本文件GeometryDef节点格式异常!', 0
      print, '样本文件GeometryDef节点格式异常!'
      return, 0
    endif
    ogeometrydef = ogeometrydefs->item(0)
    
    ocoors = ogeometrydef->GetElementsByTagName('CoordSysStr')
    if ocoors->GetLength() ne 1 then begin
      obj_destroy, [odoc, oxmlele, oregions, oregion, ogeometrydefs, ocoors]
      error = '样本文件CoordSysStr节点格式异常!'
      writelog, '样本文件CoordSysStr节点格式异常!', 0
      print, '样本文件CoordSysStr节点格式异常!'
      return, 0
    endif 
    coordinate_system = ((ocoors->item(0))->GetFirstChild())->GetNodeValue()
    
    opolygons = ogeometrydef->GetElementsByTagName('Polygon')
    num_polygon = opolygons->GetLength()
    if opolygons->GetLength() le 0 then begin
      obj_destroy, [odoc, oxmlele, oregions, oregion, ogeometrydefs, ocoors, opolygons]
      error = '样本文件Polygon节点格式异常!'
      writelog, '样本文件Polygon节点格式异常!', 0
      print, '样本文件Polygon节点格式异常!
      return, 0
    endif
    
    polygon_points = ptrarr(num_polygon)
    for j=0, num_polygon-1 do begin
      opolygon = opolygons->item(j)
      polygon_exters = opolygon->GetElementsByTagName('Exterior')
      num_exter = polygon_exters->GetLength()
      if num_exter ne 1 then begin
        obj_destroy, [odoc, oxmlele, oregions, oregion, ogeometrydefs, ocoors, opolygons, $
          opolygon, polygon_exters]
        error = 'Exterior节点格式异常!'
        writelog, 'Exterior节点格式异常!', 0
        print, 'Exterior节点格式异常!'
        return, 0
      endif
      polygon_exter = polygon_exters->item(0)
      
      polygon_linears = polygon_exter->GetElementsByTagName('LinearRing')
      num_linear = polygon_linears->GetLength()
      if num_linear ne 1 then begin
        obj_destroy, [odoc, oxmlele, oregions, oregion, ogeometrydefs, ocoors, opolygons, $
          opolygon, polygon_exters, polygon_linears]  
        error = 'LineaRing节点格式异常!'
        writelog, 'LineaRing节点格式异常!', 0
        print, 'LinearRing节点格式异常!'
        return, 0
      endif
      polygon_linear = polygon_linears->item(0)
      
      polygon_coordinates = polygon_linear->GetElementsByTagName('Coordinates')
      num_coordinates = polygon_coordinates->GetLength()
      if num_coordinates ne 1 then begin
        obj_destroy, [odoc, oxmlele, oregions, oregion, ogeometrydefs, ocoors, opolygons, $
          opolygon, polygon_exters, polygon_linears, polygon_coordinates]    
        error = 'Coordinates节点格式异常!'
        writelog, 'Coordinates节点格式异常!', 0
        print, 'Coordinates节点格式异常!'
        return, 0
      endif
      polygon_coordinate = polygon_coordinates->item(0)
      coordinates = (polygon_coordinate->GetFirstChild())->GetNodeValue()
      coordinates_ = double(call_function('strsplit', strtrim(coordinates,2), /extract))
      coordinates_num = n_elements(coordinates_)/2
      polygon_points[j] = ptr_new(reform(coordinates_, 2, coordinates_num))   
      
      obj_destroy, [opolygon,polygon_exters,polygon_exter,polygon_linears, $
                    polygon_linear,polygon_coordinates,polygon_coordinate]
    endfor    
    points[i] = ptr_new(polygon_points)  
    
    obj_destroy, [oregion, ogeometrydef, ocoors, opolygons]
  endfor
  
  obj_destroy, odoc 
  return, 1
end


function Supervised_Classification::CreateROI, $
                                 unclassifiedfile=unclassifiedfile, $
                                 points = points, $
                                 npts_per = npts_per, $
                                 error = error
  compile_opt idl2
  
  catch, error_status
  if error_status ne 0 then begin
    error = !ERROR_STATE.MSG
    WriteLog, !ERROR_STATE.MSG, 0
    catch, /cancel
    return, 0
  endif
  
  call_procedure, 'envi_open_file',unclassifiedfile,r_fid=tfid
  if tfid eq -1 then begin
    error = unclassifiedfile+'打开失败!'
    WriteLog, unclassifiedfile+'打开失败!', 0
    print, unclassifiedfile + '打开失败!'
    return, 0
  endif
  
  nRoi = n_elements(points)
  
  xfs = [0] & yfs = [0]
  for i=0, nRoi-1 do begin
    xcoordinates = reform((*points[i])[0,*])
    ycoordinates = reform((*points[i])[1,*])
    
    mapinfo = call_function('envi_get_map_info', fid=tfid, undefined=undefined)
    if undefined ne 1 then begin
      call_procedure, 'envi_convert_file_coordinates', tfid, $
        xf, yf, xcoordinates, ycoordinates      
    endif else begin
      xf = xcoordinates
      yf = ycoordinates
    endelse
    
    call_procedure, 'envi_file_query', tfid, ns=ns, nl=nl
    temp_roi_id = call_function('envi_create_roi', ns=ns, nl=nl)
    call_procedure, 'envi_define_roi', temp_roi_id, /polygon, xpts=round(xf), ypts=round(yf)
    
    addr = call_function('envi_get_roi', temp_roi_id)
    y = addr/ns
    x = addr - y*ns
    
    xfs = [xfs, x]
    yfs = [yfs, y]
    
    call_procedure, 'envi_delete_rois', temp_roi_id
  endfor
  
  roi_id = call_function('envi_create_roi', ns=ns, nl=nl)
  call_procedure, 'envi_define_roi', roi_id, /point, xpts=xfs[1:*], ypts=yfs[1:*]
  call_procedure, 'envi_get_roi_information', roi_id, npts=npts_per
  
  call_procedure, 'envi_file_mng', id=tfid, /remove
  return, roi_id  
end                           


function Supervised_Classification::StatROIInfo, $
                           unclassifiedfile = unclassifiedfile, $
                           xmlfile = xmlfile, $
                           num_class = num_class, $
                           lookup = lookup, $
                           class_names = class_names, $
                           npts = npts, $
                           covs = covs, $
                           means = means, $
                           stdvs = stdvs, $
                           error = error                        
  compile_opt idl2
  
  catch, error_status
  if error_status ne 0 then begin
    error = !ERROR_STATE.MSG
    WriteLog, !ERROR_STATE.MSG, 0
    catch, /cancel
    return, 0
  endif
  
  writelog, '开始读取ROI XML样本文件...', 1
  if ~file_test(xmlfile) then begin
    error = xmlfile+'样本文件不存在!'
    WriteLog, xmlfile+'样本文件不存在!', 0
    print, xmlfile+'样本文件不存在!'
    return, 0
  endif
  ret = self.ReadROIXML(xmlfile = xmlfile, $
                        num_class = num_class, $
                        class_names = class_names, $
                        lookup = lookup, $
                        points = points, $
                        error = error)
  if ret ne 1 then begin
    error = 'ROI XML样本文件读取失败! '+error
    writelog, 'ROI XML样本文件读取失败!', 0    
    print, '读取'+file_basename(xmlfile)+'文件失败!'
    return, 0
  endif else begin
    writelog, 'ROI XML样本文件读取完成!', 2
  endelse
  
  writelog, '样本文件开始生成ROI...', 1
  
  lookup = [[0,0,0], [lookup]]
  class_names = ['未分类', class_names]
  
  roi_ids = lonarr(num_class)
  npts = lonarr(num_class)
  for i=0, num_class-1 do begin
    res = self->CreateROI(unclassifiedfile=unclassifiedfile, points=reform(*(points[i])), npts_per=npts_per, error=error)
    if res eq 0 then begin
      error = class_names[i+1] + '类别生成ROI失败! '+error
      WriteLog, class_names[i+1] + '类别生成ROI失败!', 0
      print, class_names[i+1] + '类别生成ROI失败!'
      return, 0
    endif else begin
      roi_ids[i] = res
    endelse
    npts[i] = npts_per
  endfor
  
  writelog, '样本文件生成ROI完成!', 2
  
  writelog, '开始统计ROI相关信息...', 1
  
  call_procedure, 'envi_open_file', unclassifiedfile, r_fid=tfid
  if tfid eq -1 then begin
    error = unclassifiedfile + '文件打开失败!'
    WriteLog, unclassifiedfile + '文件打开失败!', 0
    print, unclassifiedfile + '文件打开失败!'
    return, 0
  endif
  call_procedure, 'envi_file_query', tfid, dims=dims, ns=ns, nl=nl, nb=nb
  pos = indgen(nb)  
  
  means = fltarr(nb, num_class)
  stdvs = fltarr(nb, num_class)
  covs = fltarr(nb, nb, num_class)
  
  for j=0, num_class-1 do begin
    roi_dims = [call_function('envi_get_roi_dims_ptr', roi_ids[j]), 0, 0, 0, 0]
    
    call_procedure, 'envi_doit', 'envi_stats_doit', fid=tfid, pos=pos, $
      dims=roi_dims, comp_flag=1, mean=c_mean, stdv=c_stdv
      
    call_procedure, 'envi_doit', 'envi_stats_doit', fid=tfid, pos=pos, $
      dims=roi_dims, comp_flag=4, cov=c_cov
    
    means[*,j] = c_mean
    stdvs[*,j] = c_stdv
    covs[*,*,j] = c_cov
  endfor
  
  call_procedure, 'envi_delete_rois', /all
  call_procedure, 'envi_file_mng', id=tfid, /remove
  
  writelog, 'ROI相关信息统计完成!', 2
  
  return, 1
end


function Supervised_Classification::Classify, paramfile=paramfile, error=error
  compile_opt idl2
  envi,/restore_base_save_files
  envi_batch_init
  envi_batch_status_window,/off
  
  currentdir = file_dirname(routine_filepath(/either))

  ;临时路径
  tempdir = file_dirname(currentdir)+path_sep()+'temp'+path_sep()+strtrim(string(systime(/julian),format='(f15.4)'),2)+path_sep()
  if ~file_test(tempdir) then file_mkdir, tempdir
  defsysv, '!tempdir', tempdir

  ;日志文件
  logdir = file_dirname(currentdir)+path_sep()+'log'+path_sep()
  if ~file_test(logdir) then file_mkdir, logdir
  datestr = string(systime(/julian), format='(c(CYI,":",CMOI,":",CDI))')
  strs = strsplit(datestr,':',/extract)
  year = string(strs[0], format='(I04)')
  month = string(strs[1], format='(I02)')
  day = string(strs[2], format='(I02)')
  logfile = logdir+'Supervised_Classification_'+year+month+day+'.log'
  defsysv, '!logfile', logfile

  openw, loglun, !logfile, /get_lun, /append
  printf, loglun, '**********************************************************************************'
  printf, loglun, ''
  printf, loglun, '监督分类'
  printf, loglun, 'Start time: '+string(systime(/julian),format='(c(CYI,"-",CMOI,"-",CDI," ",CHI,":",CMI,":",CSI))')
  printf, loglun, '----------------------------------------------------------'
  free_lun, loglun

  ;错误捕捉
  catch, error_status
  if error_status ne 0 then begin
    error = !ERROR_STATE.MSG
    writelog, !ERROR_STATE.MSG, 0
    catch, /cancel
    deletefiles, !tempdir
    return, 0
  endif  
  
;  if ~keyword_set(xmlfile) then xmlfile = self.xmlfile
  if ~keyword_set(paramfile) then paramfile = self.paramfile
  
  writelog, '开始读取参数据文件...', 1
  ret = self->ReadParamfile(paramfile = paramfile, $
                            xmlfile = xmlfile, $
                            method = method, $
                            unclassifiedfile = unclassifiedfile, $
                            classifiedfile = classifiedfile, $
                            error = error)  
  if ret eq 1 then begin
    writelog, file_basename(paramfile)+'参数文件读取完成!', 2
  endif else begin    
    error = file_basename(paramfile)+'参数文件读取失败! '+error
    writelog, file_basename(paramfile)+'参数文件读取失败!', 0
    deletefiles, !tempdir
;    msg = dialog_message(file_basename(paramfile)+'参数文件读取失败!', /error)
    return, 0
  endelse

  if ~file_test(unclassifiedfile) then begin
    error = unclassifiedfile+'待分类数据文件不存在!'
    writelog, unclassifiedfile+'待分类数据文件不存在!', 0
    deletefiles, !tempdir
;    msg = dialog_message(unclassifiedfile+'待分类数据文件不存在!', /error)
    return, 0
  endif
  
  writelog, '开始读取统计ROI样本信息...', 1
  ret = self->StatROIInfo(xmlfile = xmlfile, $
                          unclassifiedfile = unclassifiedfile, $
                          num_class = num_class, $
                          lookup = lookup, $
                          class_names = class_names, $
                          npts = npts, $
                          covs = covs, $
                          means = means, $
                          stdvs = stdvs, $
                          error = error)  
  if ret eq 1 then begin
    writelog, file_basename(xmlfile)+'ROI样本文件读取统计完成!', 2
  endif else begin
    error = file_basename(xmlfile)+'ROI样本文件读取统计失败! '+error
    writelog, file_basename(xmlfile)+'ROI样本文件读取统计失败!', 0
    deletefiles, !tempdir
;    msg = dialog_message(file_basename(xmlfile)+'ROI样本文件读取统计失败!', /error)
    return, 0
  endelse
  
  writelog, '开始执行监督分类...', 1
  
  call_procedure, 'envi_open_file', unclassifiedfile, r_fid=tfid
  if tfid eq -1 then begin
    error = file_basename(unclassifiedfile)+'待分类数据打开失败!'
    writelog, file_basename(unclassifiedfile)+'待分类数据打开失败!', 0
    deletefiles, !tempdir
;    msg = dialog_message(file_basename(unclassifiedfile)+'待分类数据打开失败!', /error)
    return, 0
  endif
  call_procedure, 'envi_file_query', tfid, dims=dims, ns=ns, nl=nl, nb=nb
  
  ;输出文件路径不存在则创建
  outdir = file_dirname(classifiedfile)
  if ~file_test(outdir) then file_mkdir, outdir
   
;  tempfile = outdir + path_sep() + file_basename(classifiedfile,'.tif') + '.dat' 
  tempfile = classifiedfile
  
  case method of
    'Minimum':begin
       bname = 'Min_Dist (' + file_basename(unclassifiedfile)+')'
       call_procedure, 'envi_doit', 'class_doit', $
                                     fid = tfid, $
                                     dims = dims, $
                                     pos = indgen(nb), $
                                     method = 1, $
                                     class_names = class_names, $
                                     lookup = lookup, $
                                     data_scale = data_scale, $
                                     mean = means, $
                                     stdv = stdvs, $
                                     std_mult = std_mult, $
                                     thresh = thresh, $
                                     out_bname = bname, $
                                     out_name = tempfile, $
                                     r_fid = rfid
    end
    'Maximum':begin
       bname = 'Max_Like (' + file_basename(unclassifiedfile)+')'
       call_procedure, 'envi_doit', 'class_doit', $
                                     fid = tfid, $
                                     dims = dims, $
                                     pos = indgen(nb), $
                                     method = 2, $
                                     class_names = class_names, $
                                     lookup = lookup, $
                                     mean = means, $
                                     cov = covs, $
                                     data_scale = data_scale, $
                                     stdv = stdvs, $
                                     thresh = thresh, $
                                     out_bname = bname, $
                                     out_name = tempfile, $
                                     r_fid = rfid 
                                      
    end
    'SAM':begin
       bname = 'SAM (' + file_basename(unclassifiedfile)+')'
       call_procedure, 'envi_doit', 'class_doit', $
                                     fid = tfid, $
                                     dims = dims, $
                                     pos = indgen(nb), $
                                     method = 3, $
                                     class_names = class_names, $
                                     lookup = lookup, $
                                     thresh = thresh, $
                                     mean = means, $
                                     out_bname = bname, $
                                     out_name = tempfile, $                                     
                                     r_fid = rfid
                                     
    end    
    'Mahalanobis':begin
       bname = 'Mahal_Dist (' + file_basename(unclassifiedfile)+')'
       call_procedure, 'envi_doit', 'class_doit', $
                                     fid = tfid, $
                                     dims = dims, $
                                     pos = indgen(nb), $
                                     method = 5, $
                                     class_names = class_names, $
                                     lookup = lookup, $
                                     mean = means, $
                                     cov = covs, $
                                     npts = npts, $
                                     thresh = thresh, $
                                     out_bname = bname, $
                                     out_name = tempfile, $
                                     r_fid = rfid
                                                                              
    end
    else:begin
       error = method+'监督分类方法设置有误!'
       print, method+'监督分类方法设置有误!'
       call_procedure, 'envi_file_mng', id=tfid, /remove
       writelog, method+'监督分类方法设置有误!', 0
       deletefiles, !tempdir
;       msg = dialog_message(method+'监督分类方法设置有误!', /error)
       return, 0
    endelse
  endcase
  
  call_procedure, 'envi_file_mng', id=tfid, /remove
  if rfid eq -1 then begin
    error = '监督分类失败!'
    writelog, '监督分类失败!', 0
    deletefiles, !tempdir
;    msg = dialog_message('监督分类失败!', /error)
    print, '分类失败!'
    return, 0
  endif else begin
    writelog, '监督分类完成!', 2
  endelse
  
;  writelog, '开始对分类结果进行格式转换...', 1
;  
;  call_procedure, 'envi_file_query', rfid, dims=dims, ns=ns, nl=nl, nb=nb, bnames=bnames  
;  call_procedure, 'envi_output_to_external_format', $
;     fid = rfid, $
;     dims = dims, $
;     pos = indgen(nb), $
;     out_bname = bnames, $
;     out_name = classifiedfile, $
;     /tiff  
;     
;   tfwfile = file_dirname(classifiedfile) + path_sep() + file_basename(classifiedfile,'.tif') + '.tfw'
;   if file_test(tfwfile) then file_delete, tfwfile     
;     
;   data = read_tiff(classifiedfile, geotiff=geotiff)
;   if file_test(classifiedfile) then file_delete, classifiedfile
;   write_tiff, classifiedfile, temporary(data), geotiff=geotiff;, red=lookup[0,*], green=lookup[1,*], blue=lookup[2,*]
;
;   call_procedure, 'envi_file_mng', id=rfid, /remove, /delete 
;   
;   writelog, '格式转换完成!', 2
   
   openw, loglun, !logfile, /get_lun, /append
   printf, loglun, ''
   printf, loglun, '----------------------------------------------------------'
   printf, loglun, 'End time: '+string(systime(/julian),format='(c(CYI,"-",CMOI,"-",CDI," ",CHI,":",CMI,":",CSI))')
   printf, loglun, ''
   printf, loglun, '**********************************************************************************'
   free_lun, loglun
   
   deletefiles, !tempdir
   
   return, 1
end


function Supervised_Classification::init;, $
;                      roixmlfile = roixmlfile, $
;                      paramfile = paramfile
  compile_opt idl2
  
;  self.xmlfile = roixmlfile
;  self.paramfile = paramfile
  
  return, 1
end


pro Supervised_Classification__define

  struct={Supervised_Classification, $
;                         xmlfile : '', $
                         paramfile : ''}

end


pro WriteLog, loginfo, mark
  compile_opt idl2

  openw,loglun,!logfile,/get_lun,/append
  if mark eq 0 then begin
    printf,loglun,''
    printf,loglun,string(systime(/julian),format='(c(CYI,"-",CMOI,"-",CDI," ",CHI,":",CMI,":",CSI))')
    printf,loglun,'ERROR: '+loginfo
    printf,loglun,''
    printf,loglun,'**********************************************************************************'
  endif

  if mark eq 1 then begin
    printf,loglun,''
    printf,loglun,'----------------------------------------------------------'
    printf,loglun,string(systime(/julian),format='(c(CYI,"-",CMOI,"-",CDI," ",CHI,":",CMI,":",CSI))')
    printf,loglun,'STATUS: '+loginfo
    printf,loglun,''
  endif

  if mark eq 2 then begin
    printf,loglun,''
    printf,loglun,string(systime(/julian),format='(c(CYI,"-",CMOI,"-",CDI," ",CHI,":",CMI,":",CSI))')
    printf,loglun,'STATUS: '+loginfo
    printf,loglun,'----------------------------------------------------------'
  endif

  free_lun,loglun
end


pro DeleteFiles, dir
  compile_opt idl2

  files=file_search(dir,'*',count=count)
  for i=0,count-1 do begin
    if file_test(files[i]) then file_delete,files[i],/quiet
  endfor
  file_delete,dir,/allow_nonexistent,/quiet

end