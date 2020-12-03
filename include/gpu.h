#ifndef postgpu_h
#define postgpu_h

#define CALL_GPU(KERNEL,EVENT)\
void Kernel::KERNEL() {\
  cudaDeviceSynchronize();\
  startTimer("cudaMalloc   ");\
  if( keysHost.size() > keysDevcSize ) {\
    if( keysDevcSize != 0 ) cudaFree(keysDevc);\
    cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) );\
    keysDevcSize = keysHost.size();\
  }\
  if( rangeHost.size() > rangeDevcSize ) {\
    if( rangeDevcSize != 0 ) cudaFree(rangeDevc);\
    cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) );\
    rangeDevcSize = rangeHost.size();\
  }\
  if( sourceHost.size() > sourceDevcSize ) {\
    if( sourceDevcSize != 0 ) cudaFree(sourceDevc);\
    cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(gpureal) );\
    sourceDevcSize = sourceHost.size();\
  }\
  if( targetHost.size() > targetDevcSize ) {\
    if( targetDevcSize != 0 ) cudaFree(targetDevc);\
    cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(gpureal) );\
    targetDevcSize = targetHost.size();\
  }\
  cudaDeviceSynchronize();\
  stopTimer("cudaMalloc   ");\
  startTimer("cudaMemcpy   ");\
  cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),      cudaMemcpyHostToDevice);\
  cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),     cudaMemcpyHostToDevice);\
  cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice);\
  cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice);\
  cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(gpureal));\
  cudaDeviceSynchronize();\
  stopTimer("cudaMemcpy   ");\
  cudaDeviceSynchronize();\
  startTimer(#EVENT);\
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    KERNEL##_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  cudaDeviceSynchronize();\
  stopTimer(#EVENT);\
  cudaDeviceSynchronize();\
  startTimer("cudaMemcpy   ");\
  cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(gpureal),cudaMemcpyDeviceToHost);\
  cudaDeviceSynchronize();\
  stopTimer("cudaMemcpy   ");\
}

#endif
