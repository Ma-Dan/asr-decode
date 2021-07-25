# asr-decode

## 从 [Kaldi](https://github.com/kaldi-asr/kaldi) 中裁剪的解码推理框架

## 实现
1. 不依赖OpenFST、OpenBLAS等库实现全部计算，便于学习和移植
2. 重现了基础的Viterbi解码(https://github.com/kaldi-asr/kaldi/blob/master/src/gmmbin/gmm-decode-simple.cc)

## 使用
```shell
./bin/main ./model/final.mdl ./model/HCLG.fst ./data/feat.ark
```
备注:
1. model文件来源于yesno的基础示例
2. feat.ark用下面方式计算
```
#从wave计算mfcc
kaldi/src/featbin/compute-mfcc-feats --config=conf/mfcc.conf scp:data/test_yesno/wav.scp ark:- | kaldi/src/featbin/copy-feats --compress=true ark:- ark,scp:test_yesno.ark,test_yesno.scp
#从mfcc计算cmvn
kaldi/src/featbin/compute-cmvn-stats --spk2utt=ark:data/test_yesno/spk2utt scp:data/test_yesno/feats.scp ark,scp:cmvn_test_yesno.ark,cmvn_test_yesno.scp
#应用cmvn到mfcc feature
kaldi/src/featbin/apply-cmvn --utt2spk=ark:data/test_yesno/split1/1/utt2spk scp:cmvn_test_yesno.scp scp:test_yesno.scp ark:- | kaldi/src/featbin/add-deltas ark:- ark:feat.ark
```

## Todo
1. 从音频计算MFCC等特征
2. 其他解码方式和声学模型并优化，实现[vosk-api](https://github.com/alphacep/vosk-api)的完整功能
