function Unsupervised_Classification::ReadParamfile, $
                                      paramfile = paramfile, $
                                      unclassifiedfile = unclassifiedfile, $
                                      classifiedfile = classifiedfile, $
                                      method = method, $ 
                                      params = params 
                                                                          
;                                      ; ��������
;                                      iterations = iterations, $ 
;                                      change_thresh = change_thresh, $
;                                      num_classes = num_classes, $
;                                      std_mult = std_mult, $
;                                      thresh = thresh, $
;                                      ; ISODATA
;                                      min_classes = min_classes, $                                                                            
;                                      iso_min_pixels = iso_min_pixels, $
;                                      iso_split_std = iso_split_std, $
;                                      iso_merge_dist = iso_merge_dist, $
;                                      iso_merge_pairs = iso_merge_pairs
  compile_opt idl2
  
  catch, error
  if error ne 0 then begin
    WriteLog, !ERROR_STATE.MSG, 0
    catch, /cancel
    return, 0
  endif
  
  if ~file_test(paramfile) then begin
    WriteLog, paramfile + '�����ļ�������!', 0
    print, paramfile + '�����ļ�������!'
    return, 0
  endif
  
  str = ''
  openr, lun, paramfile, /get_lun
  readf, lun, str
  temp = strsplit(str, '=', /extract)
  if n_elements(temp) ne 2 then begin
    WriteLog, paramfile + '�����ļ���������!', 0
    free_lun, lun
    return, 0
  endif
  unclassifiedfile = temp[1]
  
  readf, lun, str
  temp = strsplit(str, '=', /extract)
  if n_elements(temp) ne 2 then begin
    WriteLog, paramfile + '����ļ���������!', 0
    free_lun, lun
    return, 0
  endif
  classifiedfile = temp[1]
  
  readf, lun, str
  temp = strsplit(str, '=', /extract)
  if n_elements(temp) ne 2 then begin
    Writelog, paramfile + '���෽����������!', 0
    free_lun, lun
    return, 0
  endif
  method = temp[1]
  
  if method eq 'ISODATA' then begin
    strs = strarr(8) & vars = fltarr(8)
    readf, lun, strs
    for i=0, 7 do begin
      temp = strsplit(strs[i], '=', /extract)
      if n_elements(temp) ne 2 then begin
        WriteLog, paramfile + '�ļ���ʽ����!'
        free_lun, lun
        return, 0
      endif
      vars[i] = float(temp[1])
    endfor
  endif else if method eq 'K-Means' then begin
    strs = strarr(3) & vars = fltarr(3)
    readf, lun, strs
    for i=0, 2 do begin
      temp = strsplit(strs[i], '=', /extract)
      if n_elements(temp) ne 2 then begin
        WriteLog, paramfile + '�ļ���ʽ����!'
        free_lun, lun
        return, 0
      endif
      vars[i] = float(temp[1])
    endfor    
  endif else begin
    WriteLog, paramfile + '�����ļ��зǼල���෽����������!', 0
    free_lun, lun
    return, 0
  endelse
  
  free_lun, lun
  
  params = vars
                                
  return, 1
end


function Unsupervised_Classification::Classify, paramfile = paramfile  
  compile_opt idl2
  envi,/restore_base_save_files
  envi_batch_init
  envi_batch_status_window, /on

  currentdir = file_dirname(routine_filepath(/either))

  ;��ʱ·��
  tempdir = file_dirname(currentdir)+path_sep()+'temp'+path_sep()+strtrim(string(systime(/julian),format='(f15.4)'),2)+path_sep()
  if ~file_test(tempdir) then file_mkdir, tempdir
  defsysv, '!tempdir', tempdir

  ;��־�ļ�
  logdir = file_dirname(currentdir)+path_sep()+'log'+path_sep()
  if ~file_test(logdir) then file_mkdir, logdir
  datestr = string(systime(/julian), format='(c(CYI,":",CMOI,":",CDI))')
  strs = strsplit(datestr,':',/extract)
  year = string(strs[0], format='(I04)')
  month = string(strs[1], format='(I02)')
  day = string(strs[2], format='(I02)')
  logfile = logdir+'Unsupervised_Classification_'+year+month+day+'.log'
  defsysv, '!logfile', logfile

  openw, loglun, !logfile, /get_lun, /append
  printf, loglun, '**********************************************************************************'
  printf, loglun, ''
  printf, loglun, '�Ǽල����'
  printf, loglun, 'Start time: '+string(systime(/julian),format='(c(CYI,"-",CMOI,"-",CDI," ",CHI,":",CMI,":",CSI))')
  printf, loglun, '----------------------------------------------------------'
  free_lun, loglun

  ;����׽
  catch, error
  if error ne 0 then begin
    writelog, !ERROR_STATE.MSG, 0
    catch, /cancel
    deletefiles, !tempdir
    return, 0
  endif
  
  if ~keyword_set(paramfile) then paramfile = self.paramfile
  if ~file_test(paramfile) then begin
    writelog, paramfile+'�����ļ�������!', 0
    msg = dialog_message(paramfile+'�����ļ�������!', /error)
    return, 0
  endif
  
  writelog, '��ʼ��ȡ�����ļ�...', 1  
  ret = self->ReadParamfile(paramfile = paramfile, $
                            unclassifiedfile = unclassifiedfile, $
                            classifiedfile = classifiedfile, $
                            method = method, $
                            params = params)
  if ret eq 1 then begin
    WriteLog, file_basename(paramfile)+'��ȡ�������!', 2   
  endif else begin
    WriteLog, file_basename(paramfile)+'��ȡ����ʧ��!', 0
    return, 0
  endelse
  
  writelog, '��ʼִ�зǼල����...', 1
  
  call_procedure, 'envi_open_file', unclassifiedfile, r_fid=tfid
  if tfid eq -1 then begin
    WriteLog, unclassifiedfile + '�ļ���ʧ��!'
    return, 0
  endif
  
  call_procedure, 'envi_file_query', tfid, dims=dims, ns=ns, nl=nl, nb=nb
  
  outdir = file_dirname(classifiedfile)
  if ~file_test(outdir) then file_mkdir, outdir
  outfile = outdir + path_sep() + file_basename(classifiedfile,'.tif') + '.dat'
  
  if method eq 'ISODATA' then begin
    out_bname = 'ISODATA'+'('+file_basename(unclassifiedfile)+')'
    call_procedure,'envi_doit','class_doit', $
      fid = tfid, $
      dims = dims, $
      pos = indgen(nb), $
      out_bname = out_bname, $
      out_name = outfile, $
      method = 4, $
      min_classes = params[0], $
      num_classes = params[1], $
      iteration = params[2], $
      change_thresh = params[3], $
      iso_min_pixels = params[4], $
      iso_split_std = params[5], $
      iso_merge_dist = params[6], $
      iso_merge_pairs = params[2], $
      r_fid = rfid
          
  endif else if method eq 'K-Means' then begin
    out_bname = 'K-Means'+'('+file_basename(unclassifiedfile)+')'
    call_procedure,'envi_doit','class_doit', $
      fid = tfid, $
      dims = dims, $
      pos = indgen(nb), $
      out_bname = out_bname, $
      out_name = outfile, $
      method = 7, $
      num_classes = params[0], $
      iteration = params[1], $
      change_thresh = params[2], $
      r_fid = rfid
          
  endif else begin
    WriteLog, method+'���෽����������!', 0
    msg = dialog_message(method+'���෽����������!', /error)
    return, 0    
  endelse
  
  call_procedure, 'envi_file_mng', id=tfid, /remove
  if rfid eq -1 then begin
    WriteLog, '�Ǽල����ʧ��!', 0
    return, 0
  endif else begin
    writelog, '�Ǽල�������!', 2
  endelse
  
  writelog, '��ʼ���и�ʽת��...', 1
  
  call_procedure, 'envi_file_query', rfid, dims=dims, ns=ns, nl=nl, nb=nb, bnames=bnames, lookup=lookup  
  call_procedure, 'envi_output_to_external_format', $
     fid = rfid, $
     dims = dims, $
     pos = indgen(nb), $
     out_bname = bnames, $
     out_name = classifiedfile, $
     /tiff  
     
  data = read_tiff(classifiedfile, geotiff=geotiff)
  if file_test(classifiedfile) then file_delete, classifiedfile
  tfwfile = file_dirname(classifiedfile) + path_sep() + file_basename(classifiedfile,'.tif') + '.tfw'
  if file_test(tfwfile) then file_delete, tfwfile  

  write_tiff, classifiedfile, temporary(data), geotiff=geotiff, red=lookup[0,*], green=lookup[1,*], blue=lookup[2,*]

  call_procedure, 'envi_file_mng', id=rfid, /remove;, /delete
  
  writelog, '��ʽת�����!', 2

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


function Unsupervised_Classification::init, $
                       paramfile=paramfile
   compile_opt idl2
                       
   self.paramfile = paramfile
   
   return, 1 
end


pro Unsupervised_Classification__define
  compile_opt idl2

  struct={Unsupervised_Classification, $
    paramfile : ''}  
  
end