#define I_Width 640
#define I_Height 480
 
#define BLOCK_X  16
#define BLOCK_Y  16
// Exchange trick: Morgan McGuire, ShaderX 2008
#define s2(a,b)            { tmp = a; a = min(a,b); b = max(tmp,b); }
#define mn3(a,b,c)         s2(a,b); s2(a,c);
#define mx3(a,b,c)         s2(b,c); s2(a,c);
#define mnmx3(a,b,c)       mx3(a,b,c); s2(a,b);                               // 3 exchanges
#define mnmx4(a,b,c,d)     s2(a,b); s2(c,d); s2(a,c); s2(b,d);                // 4 exchanges
#define mnmx5(a,b,c,d,e)   s2(a,b); s2(c,d); mn3(a,c,e); mx3(b,d,e);          // 6 exchanges
#define mnmx6(a,b,c,d,e,f) s2(a,d); s2(b,e); s2(c,f); mn3(a,b,c); mx3(d,e,f); // 7 exchanges
#define SMEM(x,y)  smem[(x)+1][(y)+1]
#define IN(x,y)    d_in[(y)*nx + (x)]
__global__ void d_MedianFilter_char(int nx, int ny, char *d_out, char *d_in)
{
 int tx = threadIdx.x, ty = threadIdx.y;
 // guards: is at boundary?
 bool is_x_top = (tx == 0), is_x_bot = (tx == BLOCK_X-1);
 bool is_y_top = (ty == 0), is_y_bot = (ty == BLOCK_Y-1);
 __shared__ char smem[BLOCK_X+2][BLOCK_Y+2];
 // clear out shared memory (zero padding)
 if (is_x_top)           SMEM(tx-1, ty  ) = 0;
 else if (is_x_bot)      SMEM(tx+1, ty  ) = 0;
 if (is_y_top) {         SMEM(tx  , ty-1) = 0;
 if (is_x_top)       SMEM(tx-1, ty-1) = 0;
 else if (is_x_bot)  SMEM(tx+1, ty-1) = 0;
 } else if (is_y_bot) {  SMEM(tx  , ty+1) = 0;
 if (is_x_top)       SMEM(tx-1, ty+1) = 0;
 else if (is_x_bot)  SMEM(tx+1, ty+1) = 0;
 }
 // guards: is at boundary and still more image?
 int x = blockIdx.x * blockDim.x + tx;
 int y = blockIdx.y * blockDim.y + ty;
 is_x_top &= (x > 0); is_x_bot &= (x < nx - 1);
 is_y_top &= (y > 0); is_y_bot &= (y < ny - 1);
 // each thread pulls from image
 SMEM(tx  , ty  ) = IN(x  , y  ); // self
 if (is_x_top)           SMEM(tx-1, ty  ) = IN(x-1, y  );
 else if (is_x_bot)      SMEM(tx+1, ty  ) = IN(x+1, y  );
 if (is_y_top) {         SMEM(tx  , ty-1) = IN(x  , y-1);
 if (is_x_top)       SMEM(tx-1, ty-1) = IN(x-1, y-1);
 else if (is_x_bot)  SMEM(tx+1, ty-1) = IN(x+1, y-1);
 } else if (is_y_bot) {  SMEM(tx  , ty+1) = IN(x  , y+1);
 if (is_x_top)       SMEM(tx-1, ty+1) = IN(x-1, y+1);
 else if (is_x_bot)  SMEM(tx+1, ty+1) = IN(x+1, y+1);
 }
 __syncthreads();
 // pull top six from shared memory
 char v[6] = { SMEM(tx-1, ty-1), SMEM(tx  , ty-1), SMEM(tx+1, ty-1),
  SMEM(tx-1, ty  ), SMEM(tx  , ty  ), SMEM(tx+1, ty  ) };
 // with each pass, remove min and max values and add new value
 char tmp;
    mnmx6(v[0], v[1], v[2], v[3], v[4], v[5]);
    v[5] = SMEM(tx-1, ty+1); // add new contestant
    mnmx5(v[1], v[2], v[3], v[4], v[5]);
    v[5] = SMEM(tx  , ty+1);
    mnmx4(v[2], v[3], v[4], v[5]);
    v[5] = SMEM(tx+1, ty+1);
    mnmx3(v[3], v[4], v[5]);
 // pick the middle one
 d_out[y*nx + x] = v[4];
}
 
Ö÷»ú¶Ë£º
dim3 blocks(I_Width/BLOCK_X, I_Height/BLOCK_Y);
dim3 threads(BLOCK_X,BLOCK_Y);

d_MedianFilter_char<<<blocks,threads>>>(I_Width, I_Height, Img_Out_Dev, Img_In_Dev);