# CadQuery Playground

## 初期設定 (Ubuntu24.04.1 LTS / WSL2)

CQ-Editor を使用する前に毎回以下を実行してください。

```bash
source ./setup.shrc
```

初回は venv (仮想環境) に CQ-Editor がインストールされます。

また、`CQ-editor` が `$(pwd)/venv/bin/CQ-editor` のエイリアスになります。
エイリアスは端末を閉じるまで有効です。

## ファイルの開き方

初期設定を実施した後、以下を実行します (`path/to/script.cq.py` はスクリプトファイルのパス)。

```bash
CQ-editor path/to/script.cq.py
```
