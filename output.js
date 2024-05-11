var Potrace = function() {
    function Point(x, y) {
        this.x = x,
        this.y = y
    }
    function Bitmap(w, h) {
        this.w = w,
        this.h = h,
        this.size = w * h,
        this.arraybuffer = new ArrayBuffer(this.size),
        this.data = new Int8Array(this.arraybuffer)
    }
    function Path() {
        this.area = 0,
        this.len = 0,
        this.curve = {},
        this.pt = [],
        this.minX = 1e5,
        this.minY = 1e5,
        this.maxX = -1,
        this.maxY = -1
    }
    function Curve(n) {
        this.n = n,
        this.tag = new Array(n),
        this.c = new Array(3 * n),
        this.alphaCurve = 0,
        this.vertex = new Array(n),
        this.alpha = new Array(n),
        this.alpha0 = new Array(n),
        this.beta = new Array(n)
    }
    Point.prototype.copy = function() {
        return new Point(this.x,this.y)
    }
    ,
    Bitmap.prototype.at = function(x, y) {
        return x >= 0 && x < this.w && y >= 0 && y < this.h && 1 === this.data[this.w * y + x]
    }
    ,
    Bitmap.prototype.index = function(i) {
        var point = new Point;
        return point.y = Math.floor(i / this.w),
        point.x = i - point.y * this.w,
        point
    }
    ,
    Bitmap.prototype.flip = function(x, y) {
        this.at(x, y) ? this.data[this.w * y + x] = 0 : this.data[this.w * y + x] = 1
    }
    ,
    Bitmap.prototype.copy = function() {
        var bm = new Bitmap(this.w,this.h), i;
        for (i = 0; i < this.size; i++)
            bm.data[i] = this.data[i];
        return bm
    }
    ;
    var imgElement = document.createElement("img"), imgCanvas = document.createElement("canvas"), bm = null, maxWidth = 800, pathlist = [], callback, info = {
        isReady: !1,
        turnpolicy: "minority",
        turdsize: 5,
        optcurve: !0,
        alphamax: 1,
        opttolerance: .2
    };
    function loadImageFromFile(file) {
        info.isReady && clear(),
        imgElement.file = file;
        var reader = new FileReader, aImg;
        reader.onload = (aImg = imgElement,
        function(e) {
            aImg.src = e.target.result
        }
        ),
        reader.readAsDataURL(file)
    }
    function loadImageFromUrl(url) {
        info.isReady && clear(),
        imgElement.src = url
    }
    function loadImage(url, callback) {
        info.isReady && clear(),
        imgElement.onload = function(event) {
            loadCanvas(),
            loadBm(),
            process(callback)
        }
        ,
        imgElement.src = url
    }
    function setParameter(obj) {
        var key;
        for (key in obj)
            obj.hasOwnProperty(key) && (info[key] = obj[key])
    }
    function loadCanvas() {
        var aspectRatio = imgElement.width / imgElement.height;
        imgCanvas.width = imgElement.width,
        imgCanvas.height = imgElement.height,
        imgCanvas.width > 800 && (imgCanvas.width = 800,
        imgCanvas.height = 800 / aspectRatio),
        imgCanvas.height > 800 && (imgCanvas.height = 800,
        imgCanvas.width = 800 * aspectRatio);
        var ctx = imgCanvas.getContext("2d");
        ctx.fillStyle = "#ffffff",
        ctx.fillRect(0, 0, imgCanvas.width, imgCanvas.height),
        ctx.drawImage(imgElement, 0, 0, imgCanvas.width, imgCanvas.height)
    }
    function loadFromCanvas(canvas) {
        info.isReady && clear(),
        imgCanvas.width = canvas.width,
        imgCanvas.height = canvas.height;
        var ctx = imgCanvas.getContext("2d");
        ctx.clearRect(0, 0, imgCanvas.width, imgCanvas.height),
        ctx.drawImage(canvas, 0, 0, imgCanvas.width, imgCanvas.height),
        bm = new Bitmap(imgCanvas.width,imgCanvas.height);
        var imgdataobj = ctx.getImageData(0, 0, bm.w, bm.h), l = imgdataobj.data.length, i, j, color;
        for (i = 0,
        j = 0; i < l; i += 4,
        j++)
            imgdataobj.data[i + 3] > 0 && (color = .2126 * imgdataobj.data[i] + .7153 * imgdataobj.data[i + 1] + .0721 * imgdataobj.data[i + 2],
            bm.data[j] = color < 128 ? 1 : 0);
        info.isReady = !0
    }
    function loadBm(imgdataobj) {
        var ctx = imgCanvas.getContext("2d");
        bm = new Bitmap(imgCanvas.width,imgCanvas.height);
        var imgdataobj, l = (imgdataobj = ctx.getImageData(0, 0, bm.w, bm.h)).data.length, i, j, color;
        for (i = 0,
        j = 0; i < l; i += 4,
        j++)
            color = .2126 * imgdataobj.data[i] + .7153 * imgdataobj.data[i + 1] + .0721 * imgdataobj.data[i + 2],
            bm.data[j] = color < 128 ? 1 : 0;
        info.isReady = !0
    }
    function bmToPathlist() {
        var bm1 = bm.copy(), currentPoint = new Point(0,0), path;
        function findNext(point) {
            for (var i = bm1.w * point.y + point.x; i < bm1.size && 1 !== bm1.data[i]; )
                i++;
            return i < bm1.size && bm1.index(i)
        }
        function majority(x, y) {
            var i, a, ct;
            for (i = 2; i < 5; i++) {
                for (ct = 0,
                a = 1 - i; a <= i - 1; a++)
                    ct += bm1.at(x + a, y + i - 1) ? 1 : -1,
                    ct += bm1.at(x + i - 1, y + a - 1) ? 1 : -1,
                    ct += bm1.at(x + a - 1, y - i) ? 1 : -1,
                    ct += bm1.at(x - i, y + a) ? 1 : -1;
                if (ct > 0)
                    return 1;
                if (ct < 0)
                    return 0
            }
            return 0
        }
        function findPath(point) {
            var path = new Path, x = point.x, y = point.y, dirx = 0, diry = 1, tmp;
            for (path.sign = bm.at(point.x, point.y) ? "+" : "-"; path.pt.push(new Point(x,y)),
            x > path.maxX && (path.maxX = x),
            x < path.minX && (path.minX = x),
            y > path.maxY && (path.maxY = y),
            y < path.minY && (path.minY = y),
            path.len++,
            x += dirx,
            y += diry,
            path.area -= x * diry,
            x !== point.x || y !== point.y; ) {
                var l = bm1.at(x + (dirx + diry - 1) / 2, y + (diry - dirx - 1) / 2)
                  , r = bm1.at(x + (dirx - diry - 1) / 2, y + (diry + dirx - 1) / 2);
                r && !l ? "right" === info.turnpolicy || "black" === info.turnpolicy && "+" === path.sign || "white" === info.turnpolicy && "-" === path.sign || "majority" === info.turnpolicy && majority(x, y) || "minority" === info.turnpolicy && !majority(x, y) ? (tmp = dirx,
                dirx = -diry,
                diry = tmp) : (tmp = dirx,
                dirx = diry,
                diry = -tmp) : r ? (tmp = dirx,
                dirx = -diry,
                diry = tmp) : l || (tmp = dirx,
                dirx = diry,
                diry = -tmp)
            }
            return path
        }
        function xorPath(path) {
            var y1 = path.pt[0].y, len = path.len, x, y, maxX, minY, i, j;
            for (i = 1; i < len; i++)
                if (x = path.pt[i].x,
                (y = path.pt[i].y) !== y1) {
                    for (minY = y1 < y ? y1 : y,
                    maxX = path.maxX,
                    j = x; j < maxX; j++)
                        bm1.flip(j, minY);
                    y1 = y
                }
        }
        for (; currentPoint = findNext(currentPoint); )
            xorPath(path = findPath(currentPoint)),
            path.area > info.turdsize && pathlist.push(path)
    }
    function processPath() {
        function Quad() {
            this.data = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        }
        function Sum(x, y, xy, x2, y2) {
            this.x = x,
            this.y = y,
            this.xy = xy,
            this.x2 = x2,
            this.y2 = y2
        }
        function mod(a, n) {
            return a >= n ? a % n : a >= 0 ? a : n - 1 - (-1 - a) % n
        }
        function xprod(p1, p2) {
            return p1.x * p2.y - p1.y * p2.x
        }
        function cyclic(a, b, c) {
            return a <= c ? a <= b && b < c : a <= b || b < c
        }
        function sign(i) {
            return i > 0 ? 1 : i < 0 ? -1 : 0
        }
        function quadform(Q, w) {
            var v = new Array(3), i, j, sum;
            for (v[0] = w.x,
            v[1] = w.y,
            v[2] = 1,
            sum = 0,
            i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    sum += v[i] * Q.at(i, j) * v[j];
            return sum
        }
        function interval(lambda, a, b) {
            var res = new Point;
            return res.x = a.x + lambda * (b.x - a.x),
            res.y = a.y + lambda * (b.y - a.y),
            res
        }
        function dorth_infty(p0, p2) {
            var r = new Point;
            return r.y = sign(p2.x - p0.x),
            r.x = -sign(p2.y - p0.y),
            r
        }
        function ddenom(p0, p2) {
            var r = dorth_infty(p0, p2);
            return r.y * (p2.x - p0.x) - r.x * (p2.y - p0.y)
        }
        function dpara(p0, p1, p2) {
            var x1, y1, x2, y2;
            return x1 = p1.x - p0.x,
            y1 = p1.y - p0.y,
            x2 = p2.x - p0.x,
            x1 * (y2 = p2.y - p0.y) - x2 * y1
        }
        function cprod(p0, p1, p2, p3) {
            var x1, y1, x2, y2;
            return x1 = p1.x - p0.x,
            y1 = p1.y - p0.y,
            x2 = p3.x - p2.x,
            x1 * (y2 = p3.y - p2.y) - x2 * y1
        }
        function iprod(p0, p1, p2) {
            var x1, y1, x2, y2;
            return x1 = p1.x - p0.x,
            y1 = p1.y - p0.y,
            x1 * (x2 = p2.x - p0.x) + y1 * (y2 = p2.y - p0.y)
        }
        function iprod1(p0, p1, p2, p3) {
            var x1, y1, x2, y2;
            return x1 = p1.x - p0.x,
            y1 = p1.y - p0.y,
            x1 * (x2 = p3.x - p2.x) + y1 * (y2 = p3.y - p2.y)
        }
        function ddist(p, q) {
            return Math.sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y))
        }
        function bezier(t, p0, p1, p2, p3) {
            var s = 1 - t
              , res = new Point;
            return res.x = s * s * s * p0.x + s * s * t * 3 * p1.x + t * t * s * 3 * p2.x + t * t * t * p3.x,
            res.y = s * s * s * p0.y + s * s * t * 3 * p1.y + t * t * s * 3 * p2.y + t * t * t * p3.y,
            res
        }
        function tangent(p0, p1, p2, p3, q0, q1) {
            var A, B, C, a, b, c, d, s, r1, r2;
            return d = (b = -2 * (A = cprod(p0, p1, q0, q1)) + 2 * (B = cprod(p1, p2, q0, q1))) * b - 4 * (a = A - 2 * B + (C = cprod(p2, p3, q0, q1))) * (c = A),
            0 === a || d < 0 ? -1 : (r2 = (-b - (s = Math.sqrt(d))) / (2 * a),
            (r1 = (-b + s) / (2 * a)) >= 0 && r1 <= 1 ? r1 : r2 >= 0 && r2 <= 1 ? r2 : -1)
        }
        function calcSums(path) {
            var i, x, y;
            path.x0 = path.pt[0].x,
            path.y0 = path.pt[0].y,
            path.sums = [];
            var s = path.sums;
            for (s.push(new Sum(0,0,0,0,0)),
            i = 0; i < path.len; i++)
                x = path.pt[i].x - path.x0,
                y = path.pt[i].y - path.y0,
                s.push(new Sum(s[i].x + x,s[i].y + y,s[i].xy + x * y,s[i].x2 + x * x,s[i].y2 + y * y))
        }
        function calcLon(path) {
            var n = path.len, pt = path.pt, dir, pivk = new Array(n), nc = new Array(n), ct = new Array(4);
            path.lon = new Array(n);
            var constraint = [new Point, new Point], cur = new Point, off = new Point, dk = new Point, foundk, i, j, k1, a, b, c, d, k = 0;
            for (i = n - 1; i >= 0; i--)
                pt[i].x != pt[k].x && pt[i].y != pt[k].y && (k = i + 1),
                nc[i] = k;
            for (i = n - 1; i >= 0; i--) {
                for (ct[0] = ct[1] = ct[2] = ct[3] = 0,
                ct[dir = (3 + 3 * (pt[mod(i + 1, n)].x - pt[i].x) + (pt[mod(i + 1, n)].y - pt[i].y)) / 2]++,
                constraint[0].x = 0,
                constraint[0].y = 0,
                constraint[1].x = 0,
                constraint[1].y = 0,
                k = nc[i],
                k1 = i; ; ) {
                    if (foundk = 0,
                    ct[dir = (3 + 3 * sign(pt[k].x - pt[k1].x) + sign(pt[k].y - pt[k1].y)) / 2]++,
                    ct[0] && ct[1] && ct[2] && ct[3]) {
                        pivk[i] = k1,
                        foundk = 1;
                        break
                    }
                    if (cur.x = pt[k].x - pt[i].x,
                    cur.y = pt[k].y - pt[i].y,
                    xprod(constraint[0], cur) < 0 || xprod(constraint[1], cur) > 0)
                        break;
                    if (Math.abs(cur.x) <= 1 && Math.abs(cur.y) <= 1 || (off.x = cur.x + (cur.y >= 0 && (cur.y > 0 || cur.x < 0) ? 1 : -1),
                    off.y = cur.y + (cur.x <= 0 && (cur.x < 0 || cur.y < 0) ? 1 : -1),
                    xprod(constraint[0], off) >= 0 && (constraint[0].x = off.x,
                    constraint[0].y = off.y),
                    off.x = cur.x + (cur.y <= 0 && (cur.y < 0 || cur.x < 0) ? 1 : -1),
                    off.y = cur.y + (cur.x >= 0 && (cur.x > 0 || cur.y < 0) ? 1 : -1),
                    xprod(constraint[1], off) <= 0 && (constraint[1].x = off.x,
                    constraint[1].y = off.y)),
                    !cyclic(k = nc[k1 = k], i, k1))
                        break
                }
                0 === foundk && (dk.x = sign(pt[k].x - pt[k1].x),
                dk.y = sign(pt[k].y - pt[k1].y),
                cur.x = pt[k1].x - pt[i].x,
                cur.y = pt[k1].y - pt[i].y,
                a = xprod(constraint[0], cur),
                b = xprod(constraint[0], dk),
                c = xprod(constraint[1], cur),
                d = xprod(constraint[1], dk),
                j = 1e7,
                b < 0 && (j = Math.floor(a / -b)),
                d > 0 && (j = Math.min(j, Math.floor(-c / d))),
                pivk[i] = mod(k1 + j, n))
            }
            for (j = pivk[n - 1],
            path.lon[n - 1] = j,
            i = n - 2; i >= 0; i--)
                cyclic(i + 1, pivk[i], j) && (j = pivk[i]),
                path.lon[i] = j;
            for (i = n - 1; cyclic(mod(i + 1, n), j, path.lon[i]); i--)
                path.lon[i] = j
        }
        function bestPolygon(path) {
            function penalty3(path, i, j) {
                var n = path.len, pt = path.pt, sums = path.sums, x, y, xy, x2, y2, k, a, b, c, s, px, py, ex, ey, r = 0;
                return j >= n && (j -= n,
                r = 1),
                0 === r ? (x = sums[j + 1].x - sums[i].x,
                y = sums[j + 1].y - sums[i].y,
                x2 = sums[j + 1].x2 - sums[i].x2,
                xy = sums[j + 1].xy - sums[i].xy,
                y2 = sums[j + 1].y2 - sums[i].y2,
                k = j + 1 - i) : (x = sums[j + 1].x - sums[i].x + sums[n].x,
                y = sums[j + 1].y - sums[i].y + sums[n].y,
                x2 = sums[j + 1].x2 - sums[i].x2 + sums[n].x2,
                xy = sums[j + 1].xy - sums[i].xy + sums[n].xy,
                y2 = sums[j + 1].y2 - sums[i].y2 + sums[n].y2,
                k = j + 1 - i + n),
                px = (pt[i].x + pt[j].x) / 2 - pt[0].x,
                py = (pt[i].y + pt[j].y) / 2 - pt[0].y,
                ey = pt[j].x - pt[i].x,
                s = (ex = -(pt[j].y - pt[i].y)) * ex * (a = (x2 - 2 * x * px) / k + px * px) + 2 * ex * ey * (b = (xy - x * py - y * px) / k + px * py) + ey * ey * (c = (y2 - 2 * y * py) / k + py * py),
                Math.sqrt(s)
            }
            var i, j, m, k, n = path.len, pen = new Array(n + 1), prev = new Array(n + 1), clip0 = new Array(n), clip1 = new Array(n + 1), seg0 = new Array(n + 1), seg1 = new Array(n + 1), thispen, best, c;
            for (i = 0; i < n; i++)
                (c = mod(path.lon[mod(i - 1, n)] - 1, n)) == i && (c = mod(i + 1, n)),
                clip0[i] = c < i ? n : c;
            for (j = 1,
            i = 0; i < n; i++)
                for (; j <= clip0[i]; )
                    clip1[j] = i,
                    j++;
            for (i = 0,
            j = 0; i < n; j++)
                seg0[j] = i,
                i = clip0[i];
            for (seg0[j] = n,
            i = n,
            j = m = j; j > 0; j--)
                seg1[j] = i,
                i = clip1[i];
            for (seg1[0] = 0,
            pen[0] = 0,
            j = 1; j <= m; j++)
                for (i = seg1[j]; i <= seg0[j]; i++) {
                    for (best = -1,
                    k = seg0[j - 1]; k >= clip1[i]; k--)
                        thispen = penalty3(path, k, i) + pen[k],
                        (best < 0 || thispen < best) && (prev[i] = k,
                        best = thispen);
                    pen[i] = best
                }
            for (path.m = m,
            path.po = new Array(m),
            i = n,
            j = m - 1; i > 0; j--)
                i = prev[i],
                path.po[j] = i
        }
        function adjustVertices(path) {
            function pointslope(path, i, j, ctr, dir) {
                for (var n = path.len, sums = path.sums, x, y, x2, xy, y2, k, a, b, c, lambda2, l, r = 0; j >= n; )
                    j -= n,
                    r += 1;
                for (; i >= n; )
                    i -= n,
                    r -= 1;
                for (; j < 0; )
                    j += n,
                    r -= 1;
                for (; i < 0; )
                    i += n,
                    r += 1;
                x = sums[j + 1].x - sums[i].x + r * sums[n].x,
                y = sums[j + 1].y - sums[i].y + r * sums[n].y,
                x2 = sums[j + 1].x2 - sums[i].x2 + r * sums[n].x2,
                xy = sums[j + 1].xy - sums[i].xy + r * sums[n].xy,
                y2 = sums[j + 1].y2 - sums[i].y2 + r * sums[n].y2,
                k = j + 1 - i + r * n,
                ctr.x = x / k,
                ctr.y = y / k,
                a = (x2 - x * x / k) / k,
                b = (xy - x * y / k) / k,
                a -= lambda2 = (a + (c = (y2 - y * y / k) / k) + Math.sqrt((a - c) * (a - c) + 4 * b * b)) / 2,
                c -= lambda2,
                Math.abs(a) >= Math.abs(c) ? 0 !== (l = Math.sqrt(a * a + b * b)) && (dir.x = -b / l,
                dir.y = a / l) : 0 !== (l = Math.sqrt(c * c + b * b)) && (dir.x = -c / l,
                dir.y = b / l),
                0 === l && (dir.x = dir.y = 0)
            }
            var m = path.m, po = path.po, n = path.len, pt = path.pt, x0 = path.x0, y0 = path.y0, ctr = new Array(m), dir = new Array(m), q = new Array(m), v = new Array(3), d, i, j, k, l, s = new Point, Q, w, dx, dy, det, min, cand, xmin, ymin, z;
            for (path.curve = new Curve(m),
            i = 0; i < m; i++)
                j = po[mod(i + 1, m)],
                j = mod(j - po[i], n) + po[i],
                ctr[i] = new Point,
                dir[i] = new Point,
                pointslope(path, po[i], j, ctr[i], dir[i]);
            for (i = 0; i < m; i++)
                if (q[i] = new Quad,
                0 === (d = dir[i].x * dir[i].x + dir[i].y * dir[i].y))
                    for (j = 0; j < 3; j++)
                        for (k = 0; k < 3; k++)
                            q[i].data[3 * j + k] = 0;
                else
                    for (v[0] = dir[i].y,
                    v[1] = -dir[i].x,
                    v[2] = -v[1] * ctr[i].y - v[0] * ctr[i].x,
                    l = 0; l < 3; l++)
                        for (k = 0; k < 3; k++)
                            q[i].data[3 * l + k] = v[l] * v[k] / d;
            for (i = 0; i < m; i++) {
                for (Q = new Quad,
                w = new Point,
                s.x = pt[po[i]].x - x0,
                s.y = pt[po[i]].y - y0,
                j = mod(i - 1, m),
                l = 0; l < 3; l++)
                    for (k = 0; k < 3; k++)
                        Q.data[3 * l + k] = q[j].at(l, k) + q[i].at(l, k);
                for (; ; ) {
                    if (0 !== (det = Q.at(0, 0) * Q.at(1, 1) - Q.at(0, 1) * Q.at(1, 0))) {
                        w.x = (-Q.at(0, 2) * Q.at(1, 1) + Q.at(1, 2) * Q.at(0, 1)) / det,
                        w.y = (Q.at(0, 2) * Q.at(1, 0) - Q.at(1, 2) * Q.at(0, 0)) / det;
                        break
                    }
                    for (Q.at(0, 0) > Q.at(1, 1) ? (v[0] = -Q.at(0, 1),
                    v[1] = Q.at(0, 0)) : Q.at(1, 1) ? (v[0] = -Q.at(1, 1),
                    v[1] = Q.at(1, 0)) : (v[0] = 1,
                    v[1] = 0),
                    d = v[0] * v[0] + v[1] * v[1],
                    v[2] = -v[1] * s.y - v[0] * s.x,
                    l = 0; l < 3; l++)
                        for (k = 0; k < 3; k++)
                            Q.data[3 * l + k] += v[l] * v[k] / d
                }
                if (dx = Math.abs(w.x - s.x),
                dy = Math.abs(w.y - s.y),
                dx <= .5 && dy <= .5)
                    path.curve.vertex[i] = new Point(w.x + x0,w.y + y0);
                else {
                    if (min = quadform(Q, s),
                    xmin = s.x,
                    ymin = s.y,
                    0 !== Q.at(0, 0))
                        for (z = 0; z < 2; z++)
                            w.y = s.y - .5 + z,
                            w.x = -(Q.at(0, 1) * w.y + Q.at(0, 2)) / Q.at(0, 0),
                            dx = Math.abs(w.x - s.x),
                            cand = quadform(Q, w),
                            dx <= .5 && cand < min && (min = cand,
                            xmin = w.x,
                            ymin = w.y);
                    if (0 !== Q.at(1, 1))
                        for (z = 0; z < 2; z++)
                            w.x = s.x - .5 + z,
                            w.y = -(Q.at(1, 0) * w.x + Q.at(1, 2)) / Q.at(1, 1),
                            dy = Math.abs(w.y - s.y),
                            cand = quadform(Q, w),
                            dy <= .5 && cand < min && (min = cand,
                            xmin = w.x,
                            ymin = w.y);
                    for (l = 0; l < 2; l++)
                        for (k = 0; k < 2; k++)
                            w.x = s.x - .5 + l,
                            w.y = s.y - .5 + k,
                            (cand = quadform(Q, w)) < min && (min = cand,
                            xmin = w.x,
                            ymin = w.y);
                    path.curve.vertex[i] = new Point(xmin + x0,ymin + y0)
                }
            }
        }
        function reverse(path) {
            var curve = path.curve, m = curve.n, v = curve.vertex, i, j, tmp;
            for (i = 0,
            j = m - 1; i < j; i++,
            j--)
                tmp = v[i],
                v[i] = v[j],
                v[j] = tmp
        }
        function smooth(path) {
            var m = path.curve.n, curve = path.curve, i, j, k, dd, denom, alpha, p2, p3, p4;
            for (i = 0; i < m; i++)
                j = mod(i + 1, m),
                k = mod(i + 2, m),
                p4 = interval(.5, curve.vertex[k], curve.vertex[j]),
                0 !== (denom = ddenom(curve.vertex[i], curve.vertex[k])) ? (dd = dpara(curve.vertex[i], curve.vertex[j], curve.vertex[k]) / denom,
                alpha = (dd = Math.abs(dd)) > 1 ? 1 - 1 / dd : 0,
                alpha /= .75) : alpha = 4 / 3,
                curve.alpha0[j] = alpha,
                alpha >= info.alphamax ? (curve.tag[j] = "CORNER",
                curve.c[3 * j + 1] = curve.vertex[j],
                curve.c[3 * j + 2] = p4) : (alpha < .55 ? alpha = .55 : alpha > 1 && (alpha = 1),
                p2 = interval(.5 + .5 * alpha, curve.vertex[i], curve.vertex[j]),
                p3 = interval(.5 + .5 * alpha, curve.vertex[k], curve.vertex[j]),
                curve.tag[j] = "CURVE",
                curve.c[3 * j + 0] = p2,
                curve.c[3 * j + 1] = p3,
                curve.c[3 * j + 2] = p4),
                curve.alpha[j] = alpha,
                curve.beta[j] = .5;
            curve.alphacurve = 1
        }
        function optiCurve(path) {
            function Opti() {
                this.pen = 0,
                this.c = [new Point, new Point],
                this.t = 0,
                this.s = 0,
                this.alpha = 0
            }
            function opti_penalty(path, i, j, res, opttolerance, convc, areac) {
                var m = path.curve.n, curve = path.curve, vertex = curve.vertex, k, k1, k2, conv, i1, area, alpha, d, d1, d2, p0, p1, p2, p3, pt, A, R, A1, A2, A3, A4, s, t;
                if (i == j)
                    return 1;
                if (k = i,
                i1 = mod(i + 1, m),
                0 === (conv = convc[k1 = mod(k + 1, m)]))
                    return 1;
                for (d = ddist(vertex[i], vertex[i1]),
                k = k1; k != j; k = k1) {
                    if (k1 = mod(k + 1, m),
                    k2 = mod(k + 2, m),
                    convc[k1] != conv)
                        return 1;
                    if (sign(cprod(vertex[i], vertex[i1], vertex[k1], vertex[k2])) != conv)
                        return 1;
                    if (iprod1(vertex[i], vertex[i1], vertex[k1], vertex[k2]) < d * ddist(vertex[k1], vertex[k2]) * -.999847695156)
                        return 1
                }
                if (p0 = curve.c[3 * mod(i, m) + 2].copy(),
                p1 = vertex[mod(i + 1, m)].copy(),
                p2 = vertex[mod(j, m)].copy(),
                p3 = curve.c[3 * mod(j, m) + 2].copy(),
                area = areac[j] - areac[i],
                area -= dpara(vertex[0], curve.c[3 * i + 2], curve.c[3 * j + 2]) / 2,
                i >= j && (area += areac[m]),
                A1 = dpara(p0, p1, p2),
                A2 = dpara(p0, p1, p3),
                A3 = dpara(p0, p2, p3),
                A2 == A1)
                    return 1;
                if (s = A2 / (A2 - A1),
                0 === (A = A2 * (t = A3 / (A3 - (A4 = A1 + A3 - A2))) / 2))
                    return 1;
                for (R = area / A,
                alpha = 2 - Math.sqrt(4 - R / .3),
                res.c[0] = interval(t * alpha, p0, p1),
                res.c[1] = interval(s * alpha, p3, p2),
                res.alpha = alpha,
                res.t = t,
                res.s = s,
                p1 = res.c[0].copy(),
                p2 = res.c[1].copy(),
                res.pen = 0,
                k = mod(i + 1, m); k != j; k = k1) {
                    if (k1 = mod(k + 1, m),
                    (t = tangent(p0, p1, p2, p3, vertex[k], vertex[k1])) < -.5)
                        return 1;
                    if (pt = bezier(t, p0, p1, p2, p3),
                    0 === (d = ddist(vertex[k], vertex[k1])))
                        return 1;
                    if (d1 = dpara(vertex[k], vertex[k1], pt) / d,
                    Math.abs(d1) > opttolerance)
                        return 1;
                    if (iprod(vertex[k], vertex[k1], pt) < 0 || iprod(vertex[k1], vertex[k], pt) < 0)
                        return 1;
                    res.pen += d1 * d1
                }
                for (k = i; k != j; k = k1) {
                    if (k1 = mod(k + 1, m),
                    (t = tangent(p0, p1, p2, p3, curve.c[3 * k + 2], curve.c[3 * k1 + 2])) < -.5)
                        return 1;
                    if (pt = bezier(t, p0, p1, p2, p3),
                    0 === (d = ddist(curve.c[3 * k + 2], curve.c[3 * k1 + 2])))
                        return 1;
                    if (d1 = dpara(curve.c[3 * k + 2], curve.c[3 * k1 + 2], pt) / d,
                    d2 = dpara(curve.c[3 * k + 2], curve.c[3 * k1 + 2], vertex[k1]) / d,
                    (d2 *= .75 * curve.alpha[k1]) < 0 && (d1 = -d1,
                    d2 = -d2),
                    d1 < d2 - opttolerance)
                        return 1;
                    d1 < d2 && (res.pen += (d1 - d2) * (d1 - d2))
                }
                return 0
            }
            var curve = path.curve, m = curve.n, vert = curve.vertex, pt = new Array(m + 1), pen = new Array(m + 1), len = new Array(m + 1), opt = new Array(m + 1), om, i, j, r, o = new Opti, p0, i1, area, alpha, ocurve, s, t, convc = new Array(m), areac = new Array(m + 1);
            for (i = 0; i < m; i++)
                "CURVE" == curve.tag[i] ? convc[i] = sign(dpara(vert[mod(i - 1, m)], vert[i], vert[mod(i + 1, m)])) : convc[i] = 0;
            for (area = 0,
            areac[0] = 0,
            p0 = curve.vertex[0],
            i = 0; i < m; i++)
                i1 = mod(i + 1, m),
                "CURVE" == curve.tag[i1] && (area += .3 * (alpha = curve.alpha[i1]) * (4 - alpha) * dpara(curve.c[3 * i + 2], vert[i1], curve.c[3 * i1 + 2]) / 2,
                area += dpara(p0, curve.c[3 * i + 2], curve.c[3 * i1 + 2]) / 2),
                areac[i + 1] = area;
            for (pt[0] = -1,
            pen[0] = 0,
            len[0] = 0,
            j = 1; j <= m; j++)
                for (pt[j] = j - 1,
                pen[j] = pen[j - 1],
                len[j] = len[j - 1] + 1,
                i = j - 2; i >= 0 && !(r = opti_penalty(path, i, mod(j, m), o, info.opttolerance, convc, areac)); i--)
                    (len[j] > len[i] + 1 || len[j] == len[i] + 1 && pen[j] > pen[i] + o.pen) && (pt[j] = i,
                    pen[j] = pen[i] + o.pen,
                    len[j] = len[i] + 1,
                    opt[j] = o,
                    o = new Opti);
            for (ocurve = new Curve(om = len[m]),
            s = new Array(om),
            t = new Array(om),
            j = m,
            i = om - 1; i >= 0; i--)
                pt[j] == j - 1 ? (ocurve.tag[i] = curve.tag[mod(j, m)],
                ocurve.c[3 * i + 0] = curve.c[3 * mod(j, m) + 0],
                ocurve.c[3 * i + 1] = curve.c[3 * mod(j, m) + 1],
                ocurve.c[3 * i + 2] = curve.c[3 * mod(j, m) + 2],
                ocurve.vertex[i] = curve.vertex[mod(j, m)],
                ocurve.alpha[i] = curve.alpha[mod(j, m)],
                ocurve.alpha0[i] = curve.alpha0[mod(j, m)],
                ocurve.beta[i] = curve.beta[mod(j, m)],
                s[i] = t[i] = 1) : (ocurve.tag[i] = "CURVE",
                ocurve.c[3 * i + 0] = opt[j].c[0],
                ocurve.c[3 * i + 1] = opt[j].c[1],
                ocurve.c[3 * i + 2] = curve.c[3 * mod(j, m) + 2],
                ocurve.vertex[i] = interval(opt[j].s, curve.c[3 * mod(j, m) + 2], vert[mod(j, m)]),
                ocurve.alpha[i] = opt[j].alpha,
                ocurve.alpha0[i] = opt[j].alpha,
                s[i] = opt[j].s,
                t[i] = opt[j].t),
                j = pt[j];
            for (i = 0; i < om; i++)
                i1 = mod(i + 1, om),
                ocurve.beta[i] = s[i] / (s[i] + t[i1]);
            ocurve.alphacurve = 1,
            path.curve = ocurve
        }
        Quad.prototype.at = function(x, y) {
            return this.data[3 * x + y]
        }
        ;
        for (var i = 0; i < pathlist.length; i++) {
            var path = pathlist[i];
            calcSums(path),
            calcLon(path),
            bestPolygon(path),
            adjustVertices(path),
            "-" === path.sign && reverse(path),
            smooth(path),
            info.optcurve && optiCurve(path)
        }
    }
    function process(c) {
        c && (callback = c),
        info.isReady ? (bmToPathlist(),
        processPath(),
        callback(),
        callback = null) : setTimeout((function() {
            process(c)
        }
        ), 100)
    }
    function clear() {
        bm = null,
        pathlist = [],
        callback = null,
        info.isReady = !1
    }
    function getSVG(size, opt_type, fillcolor) {
        function path(curve) {
            function bezier(i) {
                var b = "C " + (curve.c[3 * i + 0].x * size).toFixed(3) + " " + (curve.c[3 * i + 0].y * size).toFixed(3) + ",";
                return b += (curve.c[3 * i + 1].x * size).toFixed(3) + " " + (curve.c[3 * i + 1].y * size).toFixed(3) + ",",
                b += (curve.c[3 * i + 2].x * size).toFixed(3) + " " + (curve.c[3 * i + 2].y * size).toFixed(3) + " "
            }
            function segment(i) {
                var s = "L " + (curve.c[3 * i + 1].x * size).toFixed(3) + " " + (curve.c[3 * i + 1].y * size).toFixed(3) + " ";
                return s += (curve.c[3 * i + 2].x * size).toFixed(3) + " " + (curve.c[3 * i + 2].y * size).toFixed(3) + " "
            }
            var n = curve.n, i, p = "M" + (curve.c[3 * (n - 1) + 2].x * size).toFixed(3) + " " + (curve.c[3 * (n - 1) + 2].y * size).toFixed(3) + " ";
            for (i = 0; i < n; i++)
                "CURVE" === curve.tag[i] ? p += bezier(i) : "CORNER" === curve.tag[i] && (p += segment(i));
            return p
        }
        var w = bm.w * size, h = bm.h * size, len = pathlist.length, c, i, strokec, fillc, fillrule, svg = '<svg id="svg" version="1.1" width="' + w + '" height="' + h + '" xmlns="http://www.w3.org/2000/svg">';
        for (svg += '<path d="',
        i = 0; i < len; i++)
            svg += path(c = pathlist[i].curve);
        return "curve" === opt_type ? (strokec = "black",
        fillc = "none",
        fillrule = "") : (strokec = "none",
        fillc = fillcolor,
        fillrule = ' fill-rule="evenodd"'),
        svg += '" stroke="' + strokec + '" fill="' + fillc + '"' + fillrule + "/></svg>"
    }
    function getPath(size, opt_type, fillcolor, id) {
        function path(curve) {
            function bezier(i) {
                var b = "C " + (curve.c[3 * i + 0].x * size).toFixed(3) + " " + (curve.c[3 * i + 0].y * size).toFixed(3) + ",";
                return b += (curve.c[3 * i + 1].x * size).toFixed(3) + " " + (curve.c[3 * i + 1].y * size).toFixed(3) + ",",
                b += (curve.c[3 * i + 2].x * size).toFixed(3) + " " + (curve.c[3 * i + 2].y * size).toFixed(3) + " "
            }
            function segment(i) {
                var s = "L " + (curve.c[3 * i + 1].x * size).toFixed(3) + " " + (curve.c[3 * i + 1].y * size).toFixed(3) + " ";
                return s += (curve.c[3 * i + 2].x * size).toFixed(3) + " " + (curve.c[3 * i + 2].y * size).toFixed(3) + " "
            }
            var n = curve.n, i, p = "M" + (curve.c[3 * (n - 1) + 2].x * size).toFixed(3) + " " + (curve.c[3 * (n - 1) + 2].y * size).toFixed(3) + " ";
            for (i = 0; i < n; i++)
                "CURVE" === curve.tag[i] ? p += bezier(i) : "CORNER" === curve.tag[i] && (p += segment(i));
            return p
        }
        var w = bm.w * size, h = bm.h * size, len = pathlist.length, c, i, strokec, fillc, fillrule, svgpath = "";
        for (svgpath += '<path d="',
        i = 0; i < len; i++)
            svgpath += path(c = pathlist[i].curve);
        return "curve" === opt_type ? (strokec = "black",
        fillc = "none",
        fillrule = "") : (strokec = "none",
        fillc = fillcolor,
        fillrule = ' fill-rule="evenodd"'),
        svgpath += '" stroke="' + strokec + '" fill="' + fillc + '"' + fillrule + "/>",
        id && (svgpath = svgpath.replace("<path d", '<path id="' + id + '" d')),
        svgpath
    }
    function getPathd(size, opt_type, fillcolor) {
        function path(curve) {
            function bezier(i) {
                var b = "C " + (curve.c[3 * i + 0].x * size).toFixed(3) + " " + (curve.c[3 * i + 0].y * size).toFixed(3) + ",";
                return b += (curve.c[3 * i + 1].x * size).toFixed(3) + " " + (curve.c[3 * i + 1].y * size).toFixed(3) + ",",
                b += (curve.c[3 * i + 2].x * size).toFixed(3) + " " + (curve.c[3 * i + 2].y * size).toFixed(3) + " "
            }
            function segment(i) {
                var s = "L " + (curve.c[3 * i + 1].x * size).toFixed(3) + " " + (curve.c[3 * i + 1].y * size).toFixed(3) + " ";
                return s += (curve.c[3 * i + 2].x * size).toFixed(3) + " " + (curve.c[3 * i + 2].y * size).toFixed(3) + " "
            }
            var n = curve.n, i, p = "M" + (curve.c[3 * (n - 1) + 2].x * size).toFixed(3) + " " + (curve.c[3 * (n - 1) + 2].y * size).toFixed(3) + " ";
            for (i = 0; i < n; i++)
                "CURVE" === curve.tag[i] ? p += bezier(i) : "CORNER" === curve.tag[i] && (p += segment(i));
            return p
        }
        var w = bm.w * size, h = bm.h * size, len = pathlist.length, c, i, strokec, fillc, fillrule, svgpath = "";
        for (i = 0; i < len; i++)
            svgpath += path(c = pathlist[i].curve);
        return "curve" === opt_type ? (strokec = "black",
        fillc = "none",
        fillrule = "") : (strokec = "none",
        fillc = fillcolor,
        fillrule = ' fill-rule="evenodd"'),
        svgpath
    }
    return imgElement.onload = function() {
        loadCanvas(),
        loadBm()
    }
    ,
    {
        loadImageFromFile: loadImageFromFile,
        loadImageFromUrl: loadImageFromUrl,
        loadImage: loadImage,
        loadBm: loadBm,
        setParameter: setParameter,
        loadFromCanvas: loadFromCanvas,
        process: process,
        getSVG: getSVG,
        getPath: getPath,
        getPathd: getPathd,
        img: imgElement
    }
}(Potrace);
try {
    module.exports = Potrace
} catch (err) {}
