function stem(o, R0, R, l0, l, th, sr, sf, tex_off, amb){
   let lx = .58*(o[0] + o[1] + o[2]), ly = .58*(o[3] + o[4] + o[5]),
     lz = .58*(o[6] + o[7] + o[8])
   tex_off *= tex_st
   let t = off/6, te = tex_st/sf,  dfi = 2*Math.PI/sf,
     dsi = Math.sin(dfi), dco = Math.cos(dfi),
     sith = Math.sin(th), coth = Math.cos(th),
     x2 = l*sith, y2 = l0 + l*coth,
     du = 1/(sr-1),  r = R0, dr = (R0 - R)*du
   for(let u = 0; u < 1.001; u += du ){
    let B2 = u*u,  xu = B2*x2,  yu = 2*u*(1-u)*l0 + B2*y2
    let dx = (u+u)*x2,  dy = 2*((1-u-u)*l0 + u*y2),
      no = 1/Math.sqrt(dx*dx + dy*dy),  nx = dy*no, ny = -dx*no,
      si = 0, co = 1
    for(let i = 0; i < sf+1; i++ ){
      let rsi = r*si, x = xu + rsi*nx, y = yu + rsi*ny, z = r*co
      dt[off++] = o[9]  + x*o[0] + y*o[3] + z*o[6]
      dt[off++] = o[10] + x*o[1] + y*o[4] + z*o[7]
      dt[off++] = o[11] + x*o[2] + y*o[5] + z*o[8]
      dt[off++] = amb + (1-amb)*Math.max((lx*nx + ly*ny)*si + lz*co, 0)
      dt[off++] = i/sf;  dt[off++] = tex_st*(1-u) + tex_off
      let tmp = si*dco + co*dsi;  co = co*dco - si*dsi;  si = tmp
    }
    r -= dr
   }
   for(let i = 0; i < 3; i++ ){
     o[i+9] += x2*o[i] + y2*o[i+3]
     let tmp = coth*o[i] - sith*o[i+3]
     o[i+3]  = sith*o[i] + coth*o[i+3]
     o[i] = tmp
   }
   let tf = t + sf + 1
   for(let i = 0; i < sr-1; i++ ){
     for(let j = 0; j < sf; j++ ){
       ind[di++] = t++;  ind[di++] = t;  ind[di++] = tf
       ind[di++] = tf++; ind[di++] = t;  ind[di++] = tf}
     t++;  tf++}
}
function head(o, R, h, sf, tex_off, amb){
   let lx = .58*(o[0] + o[1] + o[2]), ly = .58*(o[3] + o[4] + o[5]),
     lz = .58*(o[6] + o[7] + o[8])
   let t = off/6,  sr = 4, dr = 1/(sr-1), dfi = 2*Math.PI/sf,
     dsi = Math.sin(dfi), dco = Math.cos(dfi), R5 = .5/R
   tex_off *= tex_st
   for(let i = 0; i < sr; i++ ){
    let y = h*dr*i, r = Math.sqrt(1 - i*dr)
    let nx = (r > .000001) ? R/(2*r*h) : 0
    let ny = 1/Math.sqrt(1 + nx*nx)
    r *= R;  nx *= ny
    let si = 0, co = 1
    for(let i = 0; i < sf + 1; i++ ){
      let x = r*si,  z = r*co
      dt[off++] = o[9]  + x*o[0] + y*o[3] + z*o[6]
      dt[off++] = o[10] + x*o[1] + y*o[4] + z*o[7]
      dt[off++] = o[11] + x*o[2] + y*o[5] + z*o[8]
      dt[off++] = amb + (1-amb)*Math.max(nx*(lx*si + lz*co) + ly*ny, 0)
      dt[off++] = R5*x + .5;  dt[off++] = tex_st*(R5*z + .5) + tex_off
      let tmp = si*dco + co*dsi;  co = co*dco - si*dsi;  si = tmp
    }
   }
   let tf = t + sf + 1
   for(let i = 0; i < sr-1; i++ ){
     for(let j = 0; j < sf; j++ ){
       ind[di++] = t++;  ind[di++] = t;  ind[di++] = tf
       ind[di++] = tf++; ind[di++] = t;  ind[di++] = tf}
     t++;  tf++}
}
function cone(o, R0, R, h, sf, tex_off, amb){
   let lx = .58*(o[0] + o[1] + o[2]), ly = .58*(o[3] + o[4] + o[5]),
     lz = .58*(o[6] + o[7] + o[8])
   let t = off/6,  te = tex_st/sf, dfi = 2*Math.PI/sf,
     dsi = Math.sin(dfi), dco = Math.cos(dfi), si = 0, co = 1
   tex_off *= tex_st
   for(let i = 0; i < sf+1; i++ ){
     let x = R0*si, z = R0*co
     dt[off++] = o[9]  + x*o[0] + z*o[6]
     dt[off++] = o[10] + x*o[1] + z*o[7]
     dt[off++] = o[11] + x*o[2] + z*o[8]
     let br = amb + (1-amb)*Math.max(lx*si + lz*co, 0),
        vt = i/sf
     dt[off++] = br
      dt[off++] = vt;  dt[off++] = tex_off
     x = R*si; z = R*co
     dt[off++] = o[9]  + x*o[0] + h*o[3] + z*o[6]
     dt[off++] = o[10] + x*o[1] + h*o[4] + z*o[7]
     dt[off++] = o[11] + x*o[2] + h*o[5] + z*o[8]
     dt[off++] = br
      dt[off++] = vt;  dt[off++] = tex_off + tex_st
     let tmp = si*dco + co*dsi;  co = co*dco - si*dsi;  si = tmp
   }
   for(let j = 0; j < sf; j++ ){
     ind[di++] = t++;  ind[di++] = t++;  ind[di++] = t
     ind[di++] = t;  ind[di++] = t-1;  ind[di++] = t+1}
   for(let i = 0; i < 3; i++ ) o[i+9] += h*o[i+3]
}
function patch(o, phi, scr,scy, a, tex_off, amb){
   let lx = .58*(o[0] + o[1] + o[2]), ly = .58*(o[3] + o[4] + o[5]),
     lz = .58*(o[6] + o[7] + o[8])
   let t = off/6,  a0 = 2*(.5-a)*(1-a), a1 = 4*a*(1-a), a2 = 2*a*(a-.5)
   tex_off *= tex_st
   for(let u = 0; u < sm; u++ ){
     let u5 = 5*u,
       r = (a0*pi[0][u5] + a1*pi[1][u5] + a2*pi[2][u5])*scr,
       y = (a0*pi[0][u5+1] + a1*pi[1][u5+1] + a2*pi[2][u5+1])*scy - base,
       fv = 2*(a0*pi[0][u5+2] + a1*pi[1][u5+2] + a2*pi[2][u5+2]),
       nx = a0*pi[0][u5+3] + a1*pi[1][u5+3] + a2*pi[2][u5+3],
       ny = a0*pi[0][u5+4] + a1*pi[1][u5+4] + a2*pi[2][u5+4]
     fv *= Math.min(1 + .01/r, 2)  // tip correction
     let dfi = 2*fv/(sn - 1)
     fv += phi
     let dsi = -Math.sin(dfi),  dco = Math.cos(dfi)
     let si = Math.sin(fv),  co = Math.cos(fv)
     for(let v = 0; v < sn; v++ ){
       let x = r*si, z = r*co
       dt[off++] = o[9]  + x*o[0] + y*o[3] + z*o[6]
       dt[off++] = o[10] + x*o[1] + y*o[4] + z*o[7]
       dt[off++] = o[11] + x*o[2] + y*o[5] + z*o[8]
       dt[off++] = amb + (1-amb)*Math.max(nx*(lx*si + lz*co) + ly*ny, 0)
       dt[off++] = (v+.5)*2/sn - 1;  dt[off++] = tex_st*(u+.5)/sm + tex_off
//       dt[off++] = (u+.5)/sm;  dt[off++] = tex_st*(v+.5)/sn + tex_off
       let tmp = si*dco + co*dsi
       co = co*dco - si* dsi;  si = tmp
     }
   }
   let tn = t + sn
   for(let i = 0; i < sm-1; i++ ){
     for(let j = 0; j < sn-1; j++ ){
       ind[di++] = t++;  ind[di++] = t;  ind[di++] = tn
       ind[di++] = tn++; ind[di++] = t;  ind[di++] = tn
     }
     t++;  tn++
   }
}
function patch_tooth(o, phi, scr,scy, a, tex_off, amb){
   let lx = .58*(o[0] + o[1] + o[2]), ly = .58*(o[3] + o[4] + o[5]),
     lz = .58*(o[6] + o[7] + o[8])
   let t = off/6,  a0 = 2*(.5-a)*(1-a), a1 = 4*a*(1-a), a2 = 2*a*(a-.5)
   tex_off *= tex_st
   for(let u = 0; u < sm; u++ ){
     let u5 = 5*u,
       r = (a0*pi[0][u5] + a1*pi[1][u5] + a2*pi[2][u5])*scr,
       y = (a0*pi[0][u5+1] + a1*pi[1][u5+1] + a2*pi[2][u5+1])*scy - base,
       fv = 2*(a0*pi[0][u5+2] + a1*pi[1][u5+2] + a2*pi[2][u5+2]),
       nx = a0*pi[0][u5+3] + a1*pi[1][u5+3] + a2*pi[2][u5+3],
       ny = a0*pi[0][u5+4] + a1*pi[1][u5+4] + a2*pi[2][u5+4]
     fv *= Math.min(1 + .01/r, 2)  // tip correction
     let dfi = 2*fv/(sn - 1)
     fv += phi
     let dsi = -Math.sin(dfi),  dco = Math.cos(dfi)
     let si = Math.sin(fv),  co = Math.cos(fv)
     for(let v = 0; v < sn; v++ ){
       let x = r*si, z = r*co
       dt[off++] = o[9]  + x*o[0] + y*o[3] + z*o[6]
       dt[off++] = o[10] + x*o[1] + y*o[4] + z*o[7]
       dt[off++] = o[11] + x*o[2] + y*o[5] + z*o[8]
       dt[off++] = amb + (1-amb)*Math.max(nx*(lx*si + lz*co) + ly*ny, 0)
       dt[off++] = (v+.5)*2/sn - 1;  dt[off++] = tex_st*(u+.5)/sm + tex_off
       let tmp = si*dco + co*dsi
       co = co*dco - si* dsi;  si = tmp
     }
   }
   let tn = t + sn
   for(let i = 0; i < sm-1; i++ ){
     ind[di++] = t++;  ind[di++] = t;  ind[di++] = ++tn
     for(let j = 0; j < sn-3; j++ ){
       ind[di++] = t++;  ind[di++] = t;  ind[di++] = tn
       ind[di++] = tn++; ind[di++] = t;  ind[di++] = tn
     }
       ind[di++] = t++;  ind[di++] = t;  ind[di++] = tn++
     t++;  tn++
   }
}
function patch_tooth2(o, phi, scr,scy, a, tex_off, amb){
   let lx = .58*(o[0] + o[1] + o[2]), ly = .58*(o[3] + o[4] + o[5]),
     lz = .58*(o[6] + o[7] + o[8])
   let t = off/6,  a0 = 2*(.5-a)*(1-a), a1 = 4*a*(1-a), a2 = 2*a*(a-.5)
   tex_off *= tex_st
   for(let u = 0; u < sm; u++ ){
     let u5 = 5*u,
       r = (a0*pi[0][u5] + a1*pi[1][u5] + a2*pi[2][u5])*scr,
       y = (a0*pi[0][u5+1] + a1*pi[1][u5+1] + a2*pi[2][u5+1])*scy - base,
       fv = 2*(a0*pi[0][u5+2] + a1*pi[1][u5+2] + a2*pi[2][u5+2]),
       nx = a0*pi[0][u5+3] + a1*pi[1][u5+3] + a2*pi[2][u5+3],
       ny = a0*pi[0][u5+4] + a1*pi[1][u5+4] + a2*pi[2][u5+4]
     fv *= Math.min(1 + .01/r, 2)  // tip correction
     let dfi = 2*fv/(sn - 1)
     fv += phi
     let dsi = -Math.sin(dfi),  dco = Math.cos(dfi)
     let si = Math.sin(fv),  co = Math.cos(fv)
     for(let v = 0; v < sn; v++ ){
       let x = r*si, z = r*co
       dt[off++] = o[9]  + x*o[0] + y*o[3] + z*o[6]
       dt[off++] = o[10] + x*o[1] + y*o[4] + z*o[7]
       dt[off++] = o[11] + x*o[2] + y*o[5] + z*o[8]
       dt[off++] = amb + (1-amb)*Math.max(nx*(lx*si + lz*co) + ly*ny, 0)
       dt[off++] = (v+.5)*2/sn - 1;  dt[off++] = tex_st*(u+.5)/sm + tex_off
       let tmp = si*dco + co*dsi
       co = co*dco - si* dsi;  si = tmp
     }
   }
   let tn = t + sn
   for(let j = 0; j < (sn-1)/2; j++ ){
     ind[di++] = tn++;  ind[di++] = t++;  ind[di++] = tn
     ind[di++] = tn++; ind[di++] = ++t;  ind[di++] = tn
   }
   t++;  tn++
   for(let i = 1; i < sm-1; i++ ){
     for(let j = 0; j < sn-1; j++ ){
       ind[di++] = t++;  ind[di++] = t;  ind[di++] = tn
       ind[di++] = tn++; ind[di++] = t;  ind[di++] = tn
     }
     t++;  tn++
   }
}
function bezier1DD(CP, pi){
   let dt = 1/(sm-1),  t = 0, j = 0
   for (let i = 0; i < sm; i++){
     let t2 = t*t,  t1 = 1-t,  t12 = t1*t1,  tt1 = 2*t*t1
     let b3 = t12*t1,  b2 = 3*t*t12,  b1 = 3*t2*t1,  b0 = t*t2
     let d3 = -t12,  d2 = t12 - 2*tt1,  d1 = 2*tt1 - t2,  d0 = t2
     pi[j++] = b0*CP[0] + b1*CP[3] + b2*CP[6] + b3*CP[9]
     pi[j++] = b0*CP[1] + b1*CP[4] + b2*CP[7] + b3*CP[10]
     pi[j++] = b0*CP[2] + b1*CP[5] + b2*CP[8] + b3*CP[11]
     let dr = d0*CP[0] + d1*CP[3] + d2*CP[6] + d3*CP[9],
         dy = d0*CP[1] + d1*CP[4] + d2*CP[7] + d3*CP[10],
         no = Math.sqrt(dr*dr + dy*dy)
     if(no > .000001){ pi[j++] = dy/no;  pi[j++] = -dr/no}
     else{ pi[j++] = 0;  pi[j++] = 1}
     t += dt
   }
}
function transp(s, id, tex, img){
   let t = 4*id*s*s,  s4 = 4*s
   for(let i = 0; i < s4; i += 4 ){
     let js = 0
     for(let j = 0; j < s; j++ ){
       img[t++] = tex[js + i];  img[t++] = tex[js + i + 1]
       img[t++] = tex[js + i + 2];  img[t++] = 255
       js += s4 }}
}
