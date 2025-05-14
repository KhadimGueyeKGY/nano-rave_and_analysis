#!/bin/sh

echo "🛠️ Installing tools from requirements.txt..."

while read tool; do
  if [ ! -z "$tool" ]; then
    echo "➡️ Installing $tool"
    sudo apt-get install -y "$tool"
  fi
done < requirements.txt

echo "✅ Done!"
