# Papers-Reading

## LLM

### Survey

|Date|Paper|Key Words|Github|
|:---:|:---:|:---:|:---:|
|2024.4.22|[A Survey on Efficient Inference for Large Language Models](https://arxiv.org/abs/2404.14294)|Efficient Inference||
|2024.12.27|[A Survey on Large Language Model Acceleration based on KV Cache Management](https://arxiv.org/abs/2412.19442)|KV Cache Management|[Awesome-KV-Cache-Management](https://github.com/TreeAI-Lab/Awesome-KV-Cache-Management) & [Awesome-LLM-KV-Cache](https://github.com/Zefan-Cai/Awesome-LLM-KV-Cache)|

### Models

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2023.4.17|[Visual Instruction Tuning](https://arxiv.org/abs/2304.08485)|LLaVa|
|2024.3.8|[DeepSeek-VL: Towards Real-World Vision-Language Understanding](https://arxiv.org/abs/2403.05525)|DeepSeek-VL: Dense & VLM|
|2024.7.10|[PaliGemma: A versatile 3B VLM for transfer](https://arxiv.org/abs/2407.07726)|Google small VLM: Paligemma|
|2024.10.8|[Aria: An Open Multimodal Native Mixture-of-Experts Model](https://arxiv.org/abs/2410.05993)|First MoE VLM: Aria|
|2024.12.6|[Expanding Performance Boundaries of Open-Source Multimodal Models with Model, Data, and Test-Time Scaling](https://arxiv.org/abs/2412.05271)|VLM: InternVL 2.5|
|2024.12.13|[DeepSeek-VL2: Mixture-of-Experts Vision-Language Models for Advanced Multimodal Understanding](https://arxiv.org/abs/2412.10302)|DeepSeek-VL2: MOE & VLM|
|2024.12.27|[DeepSeek-V3 Technical Report](https://arxiv.org/abs/2412.19437)|DeepSeek-V3 Technical Report|

### Kernel Optimization

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2018.11.19|[Modeling Deep Learning Accelerator Enabled GPUs](https://arxiv.org/abs/1811.08309)|Tensor Core Design && GPGPU-Sim|
|2019|[Understanding the Overheads of Launching CUDA Kernels](https://www.hpcs.cs.tsukuba.ac.jp/icpp2019/data/posters/Poster17-abst.pdf)|Kernel Launch Overhead|
|2021.10.25|[Bolt: Bridging the Gap between Auto-tuners and Hardware-native Performance](https://arxiv.org/abs/2110.15238)|Persistent kernel fusion|
|2022.4.5|[PERKS: a Locality-Optimized Execution Model for Iterative Memory-bound GPU Applications](https://arxiv.org/abs/2204.02064)|PERsistent KernelS (PERKS)|
|2023.12.19|[A Case Study in CUDA Kernel Fusion: Implementing FlashAttention-2 on NVIDIA Hopper Architecture using the CUTLASS Library](https://arxiv.org/abs/2312.11918)|FlashAttention2 using cutlass|
|2025.4.8|[Accelerating LLM Inference Throughput via Asynchronous KV Cache Prefetching](https://arxiv.org/abs/2504.06319)|Prefetches required KV Cache into GPU L2 cache|

### Serving

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2024.5.7|[QServe: W4A8KV4 Quantization and System Co-design for Efficient LLM Serving](https://arxiv.org/abs/2405.04532v2)|Boosts efficiency with W4A8KV4 quantization & Reduces dequantization overheads|
|2025.2.20|[LServe: Efficient Long-sequence LLM Serving with Unified Sparse Attention](https://arxiv.org/abs/2502.14866)|Accelerates long-context LLM inference through unified sparse attention & Hierarchical KV cache management|

### Training

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2021.7.14|[Chimera: Efficiently Training Large-Scale Neural Networks with Bidirectional Pipelines](https://arxiv.org/abs/2107.06925)|Bidirectional Pipelines|
|2023.11.30|[Zero Bubble Pipeline Parallelism](https://arxiv.org/abs/2401.10241)|Zero Bubble PP|

### Attention

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2022.5.27|[FlashAttention: Fast and Memory-Efficient Exact Attention with IO-Awareness](https://arxiv.org/abs/2205.14135)|Flash Attention|
|2023.7.18|[FlashAttention-2: Faster Attention with Better Parallelism and Work Partitioning](https://tridao.me/publications/flash2/flash2.pdf)| Flash Attention 2|
|2024.7.12|[FlashAttention-3 is optimized for Hopper GPUs (e.g. H100)](https://tridao.me/publications/flash3/flash3.pdf)|Flash Attention 3|
|2024.10.3|[SageAttention: Accurate 8-Bit Attention for Plug-and-play Inference Acceleration](https://arxiv.org/abs/2410.02367)|Sage Attention|
|2024.11.17|[SageAttention2: Efficient Attention with Thorough Outlier Smoothing and Per-thread INT4 Quantization](https://arxiv.org/abs/2411.10958)|Sage Attention 2|
|2024.3.7|[Slim attention: cut your context memory in half without loss of accuracy -- K-cache is all you need for MHA](https://arxiv.org/abs/2503.05840)|slim Attention|
|2025.4.1|[Multi-Token Attention](https://arxiv.org/abs/2504.00927v1)|Multi-Token Attention|

### Quantization

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2022.6.4|[ZeroQuant: Efficient and Affordable Post-Training Quantization for Large-Scale Transformers](https://arxiv.org/abs/2206.01861)|INT8 weights and INT8 activations|
|2022.8.15|[LLM.int8(): 8-bit Matrix Multiplication for Transformers at Scale](https://arxiv.org/abs/2208.07339)|LLM.int8|
|2022.11.18|[SmoothQuant: Accurate and Efficient Post-Training Quantization for Large Language Models](https://arxiv.org/abs/2211.10438)|8-bit Weight，8-bit Activation (W8A8)|
|2023.5.23|[Memory-Efficient Fine-Tuning of Compressed Large Language Models via sub-4-bit Integer Quantization](https://arxiv.org/abs/2305.14152)|Parameter-Efficient and Quantization-aware Adaptation (PEQA) [LLM-QAT]|
|2023.5.23|[QLoRA: Efficient Finetuning of Quantized LLMs](https://arxiv.org/abs/2305.14314)|QLoRA & NF4 (4-bit NormalFloat) [LLM-QAT]|
|2023.5.29|[LLM-QAT: Data-Free Quantization Aware Training for Large Language Models](https://arxiv.org/abs/2305.17888)|LLM Quantization Aware Training [LLM-QAT]|
|2023.3.13|[FlexGen: High-Throughput Generative Inference of Large Language Models with a Single GPU](https://arxiv.org/abs/2303.06865)|KV Cache 4-bit|
|2023.6.1|[AWQ: Activation-aware Weight Quantization for LLM Compression and Acceleration](https://arxiv.org/abs/2306.00978)|Activation-aware Weight Quantization (AWQ)|
|2023.6.13|[SqueezeLLM: Dense-and-Sparse Quantization](https://arxiv.org/abs/2306.07629)|KV Cache 3-bit|
|2024.1.31|[KVQuant: Towards 10 Million Context Length LLM Inference with KV Cache Quantization](https://arxiv.org/abs/2401.18079)|KV Cache 2、3、4-bit|
|2024.2.5|[KIVI: A Tuning-Free Asymmetric 2bit Quantization for KV Cache](https://arxiv.org/abs/2402.02750)|KV Cache 2-bit|
|2024.2.26|[A Comprehensive Evaluation of Quantization Strategies for Large Language Models](https://arxiv.org/abs/2402.16775)|PTQ|
|2024.3.8|[GEAR: An Efficient KV Cache Compression Recipe for Near-Lossless Generative Inference of LLM](https://arxiv.org/abs/2403.05527)|KV Cache Compression|
|2024.6.5|[QJL: 1-Bit Quantized JL Transform for KV Cache Quantization with Zero Overhead](https://arxiv.org/abs/2406.03482)|3 Bits KV Cache|
|2024.11.26|[Efficient LLM Inference with I/O-Aware Partial KV Cache Recomputation](https://arxiv.org/abs/2411.17089)|KV Cache Recomputation|
|2025.1.25|[RotateKV: Accurate and Robust 2-Bit KV Cache Quantization for LLMs via Outlier-Aware Adaptive Rotations](https://arxiv.org/abs/2501.16383)|2-Bit KV Cache|
|2025.2.4|[ParetoQ: Scaling Laws in Extremely Low-bit LLM Quantization](https://arxiv.org/abs/2502.02631)|Low-bit LLM Quantization|
|2025.2.15|[CalibQuant: 1-Bit KV Cache Quantization for Multimodal LLMs](https://arxiv.org/abs/2502.14882)|1-Bit KV Cache|
|2025.3.25|[LogQuant: Log-Distributed 2-Bit Quantization of KV Cache with Superior Accuracy Preservation](https://arxiv.org/abs/2503.19950)|2-Bit KV Cache|

### MOE

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2021.1.11|[Switch Transformers: Scaling to Trillion Parameter Models with Simple and Efficient Sparsity](https://arxiv.org/abs/2101.03961)|Mixture of Expert (MoE)|
|2024.1.11|[DeepSeekMoE: Towards Ultimate Expert Specialization in Mixture-of-Experts Language Models](https://arxiv.org/abs/2401.06066)|DeepSeekMoE|

### Inference

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2017.6.12|[Attention Is All You Need](https://arxiv.org/abs/1706.03762)|Transformer & Attention|
|2018.6.11|[Improving Language Understanding by Generative Pre-Training](https://cdn.openai.com/research-covers/language-unsupervised/language_understanding_paper.pdf)|Generative transformer model|
|2018.10.11|[BERT: Pre-training of Deep Bidirectional Transformers for Language Understanding](https://arxiv.org/abs/1810.04805)|BERT (Bidirectional Encoder Representations from Transformers)|
|2019.1.9|[Transformer-XL: Attentive Language Models Beyond a Fixed-Length Context](https://arxiv.org/abs/1901.02860)|Transformer-XL (extra-long)|
|2019.5.17|[ERNIE: Enhanced Language Representation with Informative Entities](https://arxiv.org/abs/1905.07129)|Knowledge graphs with BERT|
|2024.2.27|[Actions Speak Louder than Words: Trillion-Parameter Sequential Transducers for Generative Recommendations](https://arxiv.org/abs/2402.17152)|LLM for Large-scale recommendation systems|
|2024.3.19|[When Do We Not Need Larger Vision Models?](https://arxiv.org/abs/2403.13043)|Scaling on Scales|
|2024.7.28|[Enhancing Taobao Display Advertising with Multimodal Representations: Challenges, Approaches and Insights](https://arxiv.org/abs/2407.19467)|Advertising with Multimodal|
|2024.8.22|[NanoFlow: Towards Optimal Large Language Model Serving Throughput](https://arxiv.org/abs/2408.12757)|A novel serving framework: NanoFlow|

### Transformer

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2020.10.22|[An Image is Worth 16x16 Words: Transformers for Image Recognition at Scale](https://arxiv.org/abs/2010.11929)|Vision Transformer (ViT)|

### Others

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|2019.2.24|[Language Models are Unsupervised Multitask Learners](https://cdn.openai.com/better-language-models/language_models_are_unsupervised_multitask_learners.pdf)|GPT-2|
|2019.10.2|[DistilBERT, a distilled version of BERT: smaller, faster, cheaper and lighter](https://arxiv.org/abs/1910.01108)|Bert distilled version & knowledge distillation|
|2019.10.23|[Exploring the Limits of Transfer Learning with a Unified Text-to-Text Transformer](https://arxiv.org/abs/1910.10683)| Unified Text-to-Text Transformer & T5 (Encoder-Decoder)|
|2020.05.22|[Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks](https://arxiv.org/abs/2005.11401)| Retrieval-Augmented Generation (RAG)|
|2020.5.28|[Language Models are Few-Shot Learners](https://arxiv.org/abs/2005.14165)|GPT-3|
|2020.10.29|[AutoPrompt: Eliciting Knowledge from Language Models with Automatically Generated Prompts](https://arxiv.org/abs/2010.15980)|Auto generate prompt|
|2021.4.20|[RoFormer: Enhanced Transformer with Rotary Position Embedding](https://arxiv.org/abs/2104.09864)|RoPE|
|2021.6.17|[LoRA: Low-Rank Adaptation of Large Language Models](https://arxiv.org/abs/2106.09685)|LoRA|
|2021.7.7|[Evaluating Large Language Models Trained on Code](https://arxiv.org/abs/2107.03374)|finetune|
|2021.9.3|[Finetuned Language Models Are Zero-Shot Learners](https://arxiv.org/abs/2109.01652)|finetune|
|2021.12.13|[GLaM: Efficient Scaling of Language Models with Mixture-of-Experts](https://arxiv.org/abs/2112.06905)|GLaM & MOE|
|2021.12.17|[WebGPT: Browser-assisted question-answering with human feedback](https://arxiv.org/abs/2112.09332)|WebGPT|
|2025.5.14|[Insights into DeepSeek-V3: Scaling Challenges and Reflections on Hardware for AI Architectures](https://arxiv.org/abs/2505.09343v1)|DeepSeek's AI Architectures|

## Algorithm

|Date|Paper|Key Words|
|:---:|:---:|:---:|
|1972|[Karp's 21 NP-complete problems](https://en.wikipedia.org/wiki/Karp%27s_21_NP-complete_problems)|Karp's 21 NP-complete problems|
|1973|[An n^{5/2} algorithm for maximum matchings in bipartite graphs](https://web.eecs.umich.edu/~pettie/matching/Hopcroft-Karp-bipartite-matching.pdf)|Hopcroft-Karp Algorithm|
|2002|[A 27/26-Approximation Algorithm for the Chromatic Sum Coloring of Bipartite Graphs](https://eti.pg.edu.pl/documents/174618/23783336/A%202726-Approximation%20Algorithm%20for%20the%20Chromatic%20Sum%20Coloring%20of%20Bipartite%20Graphs.pdf)|Chromatic Sum Coloring of Bipartite Graphs|
|2015.6.16|[An Efficient Data Structure for Processing Palindromes in Strings](https://arxiv.org/abs/1506.04862)|Palindromic Tree|
|2017.8.11|[An Introduction to Quantum Computing, Without the Physics](https://arxiv.org/abs/1708.03684)| Quantum Computing, Without the Physics|
|2018.7.30|[A Simple Near-Linear Pseudopolynomial Time Randomized Algorithm for Subset Sum](https://arxiv.org/abs/1807.11597)|A Simple Near-Linear Pseudopolynomial Time Randomized Algorithm for Subset Sum|
|2021.2.11|[Hybrid Neural Fusion for Full-frame Video Stabilization](https://arxiv.org/abs/2102.06205)|Video Stabilization Algorithm|
|2022.11.21|[The Berlekamp-Massey Algorithm revisited](http://hlombardi.free.fr/publis/BMAvar.pdf)|Berlekamp-Massey Algorithm|

