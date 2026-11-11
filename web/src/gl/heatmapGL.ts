// WebGL2 heatmap. The cells × bins matrix is one R32F texture (rows in dataset
// order). Two small lookup textures map the screen to it:
//   - a row map: display row → dataset row (the clustering order, after filters)
//   - a column map: screen pixel column → bin (bins take their genomic width;
//     -1 draws the background, e.g. past the genome's end)
// and a colour lookup texture chosen by the layer's semantics colours the value.
// One draw call per frame, whatever the size.

import { LINEAR_LUT, type Scale } from "../lib/scales";

const VS = `#version 300 es
in vec2 aPos;
out vec2 vUv;
void main() {
  vUv = vec2((aPos.x + 1.0) * 0.5, (1.0 - aPos.y) * 0.5); // y = 0 at the top
  gl_Position = vec4(aPos, 0.0, 1.0);
}`;

const MAP_W = 4096;

const FS = `#version 300 es
precision highp float;
precision highp int;
uniform highp sampler2D uData;
uniform highp sampler2D uLut;
uniform highp sampler2D uRows;   // display row -> dataset row, MAP_W wide
uniform highp sampler2D uCols;   // pixel column -> bin (-1: none)
uniform vec2 uRowRange;          // first display row, last (exclusive), may be fractional
uniform int uNCols;
uniform int uDiscrete;
uniform float uLo, uHi;
uniform int uLutN;
uniform vec4 uNa;
uniform vec4 uBg;
in vec2 vUv;
out vec4 o;
void main() {
  int px = min(uNCols - 1, int(vUv.x * float(uNCols)));
  float bin = texelFetch(uCols, ivec2(px, 0), 0).r;
  if (bin < 0.0) { o = uBg; return; }
  int drow = int(floor(uRowRange.x + vUv.y * (uRowRange.y - uRowRange.x)));
  float row = texelFetch(uRows, ivec2(drow % ${MAP_W}, drow / ${MAP_W}), 0).r;
  float v = texelFetch(uData, ivec2(int(bin), int(row)), 0).r;
  if (v < -1e29) { o = uNa; return; }   // missing (NaN is sent as -1e30)
  int k;
  if (uDiscrete == 1) k = int(clamp(floor(v + 0.5), uLo, uHi) - uLo);
  else k = int(floor(clamp((v - uLo) / (uHi - uLo), 0.0, 1.0) * float(uLutN - 1) + 0.5));
  o = texelFetch(uLut, ivec2(k, 0), 0);
}`;

const NA = -1e30;

function compile(gl: WebGL2RenderingContext, type: number, src: string) {
  const s = gl.createShader(type)!;
  gl.shaderSource(s, src);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s) || "shader");
  return s;
}

type RGB = [number, number, number];

export class HeatmapGL {
  private gl: WebGL2RenderingContext;
  private prog: WebGLProgram;
  private tex: Record<"data" | "lut" | "rows" | "cols", WebGLTexture>;
  private u: Record<string, WebGLUniformLocation | null> = {};
  private scale: Scale | null = null;
  private lutN = 1;
  private nCols = 1;
  private dataKey: unknown = null;
  rows = 0;
  cols = 0;
  readonly maxTexture: number;

  constructor(readonly canvas: HTMLCanvasElement) {
    const gl = canvas.getContext("webgl2", { antialias: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error("WebGL2 is not available in this browser");
    this.gl = gl;
    this.maxTexture = gl.getParameter(gl.MAX_TEXTURE_SIZE);
    const p = gl.createProgram()!;
    gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, VS));
    gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, FS));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p) || "link");
    this.prog = p;
    for (const n of ["uData", "uLut", "uRows", "uCols", "uRowRange", "uNCols", "uDiscrete", "uLo", "uHi", "uLutN", "uNa", "uBg"])
      this.u[n] = gl.getUniformLocation(p, n);

    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(p, "aPos");
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);

    const mk = () => {
      const t = gl.createTexture()!;
      gl.bindTexture(gl.TEXTURE_2D, t);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      return t;
    };
    this.tex = { data: mk(), lut: mk(), rows: mk(), cols: mk() };
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
  }

  private upload(t: WebGLTexture, w: number, h: number, data: Float32Array) {
    const gl = this.gl;
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, w, h, 0, gl.RED, gl.FLOAT, data);
  }

  /** Upload a row-major rows × cols matrix (NaN = missing). Skipped if `values` is unchanged. */
  setData(values: Float32Array, rows: number, cols: number) {
    if (values === this.dataKey) return;
    if (rows > this.maxTexture || cols > this.maxTexture)
      throw new Error(`matrix ${rows}×${cols} exceeds this GPU's texture limit (${this.maxTexture})`);
    const v = new Float32Array(values.length);
    for (let i = 0; i < values.length; i++) v[i] = Number.isNaN(values[i]) ? NA : values[i];
    this.upload(this.tex.data, cols, rows, v);
    this.rows = rows;
    this.cols = cols;
    this.dataKey = values;
  }

  /** Dataset rows in display order. */
  setRows(order: ArrayLike<number>) {
    const n = Math.max(1, order.length), h = Math.ceil(n / MAP_W);
    const m = new Float32Array(MAP_W * h);
    for (let i = 0; i < order.length; i++) m[i] = order[i];
    this.upload(this.tex.rows, MAP_W, h, m);
  }

  /** Bin of each pixel column (-1 for none). */
  setCols(bins: Float32Array) {
    this.nCols = bins.length;
    this.upload(this.tex.cols, bins.length, 1, bins);
  }

  setScale(scale: Scale, lut: Uint8Array) {
    const gl = this.gl;
    this.scale = scale;
    this.lutN = scale.kind === "discrete" ? lut.length / 4 : LINEAR_LUT;
    gl.bindTexture(gl.TEXTURE_2D, this.tex.lut);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, this.lutN, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, lut);
  }

  /** Draw display rows [r0, r1) over the whole canvas. */
  draw(r0: number, r1: number, na: RGB, bg: RGB) {
    const gl = this.gl, s = this.scale;
    if (!s || !this.rows) return;
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
    gl.useProgram(this.prog);
    const bind = (unit: number, t: WebGLTexture, name: string) => {
      gl.activeTexture(gl.TEXTURE0 + unit);
      gl.bindTexture(gl.TEXTURE_2D, t);
      gl.uniform1i(this.u[name], unit);
    };
    bind(0, this.tex.data, "uData");
    bind(1, this.tex.lut, "uLut");
    bind(2, this.tex.rows, "uRows");
    bind(3, this.tex.cols, "uCols");
    gl.uniform2f(this.u.uRowRange, r0, r1);
    gl.uniform1i(this.u.uNCols, this.nCols);
    gl.uniform1i(this.u.uDiscrete, s.kind === "discrete" ? 1 : 0);
    gl.uniform1f(this.u.uLo, s.lo);
    gl.uniform1f(this.u.uHi, s.hi);
    gl.uniform1i(this.u.uLutN, this.lutN);
    gl.uniform4f(this.u.uNa, na[0] / 255, na[1] / 255, na[2] / 255, 1);
    gl.uniform4f(this.u.uBg, bg[0] / 255, bg[1] / 255, bg[2] / 255, 1);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }
}
