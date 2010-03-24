function my_write_ascii_2D, data_array, name_as_string
openw, 1, name_as_string
; Annahme: array[n,m] = array [Spalten, Zeilen], so soll auch geschrieben werden, = fltarr(m,n)
Anzahl_Spalten= n_elements(data_array[*,0])
form='('+string(Anzahl_Spalten)+'(E))'
printf, 1, FORMAT=form,data_array
close, 1
end